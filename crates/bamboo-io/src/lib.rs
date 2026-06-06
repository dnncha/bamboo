//! Unified local and cloud object I/O for Bamboo.
//!
//! Supports plain filesystem paths plus `file://`, `s3://`, `gs://`, and `https://`
//! locations via the [`object_store`] crate.

mod error;
mod runtime;

pub use error::IoError;

/// How a BAM is accessed by Bamboo readers.
#[derive(Debug, Clone)]
pub enum BamStorage {
    /// Stream from a local filesystem path (no full-file buffering on open).
    Local(std::path::PathBuf),
    /// Fully buffered object for remote URIs.
    Remote {
        uri: String,
        data: std::sync::Arc<[u8]>,
        index_data: Option<std::sync::Arc<[u8]>>,
        index_uri: Option<String>,
    },
}

/// A resolved BAM object and optional sidecar index bytes.
#[derive(Debug, Clone)]
pub struct BamSource {
    pub uri: String,
    pub storage: BamStorage,
}

impl BamSource {
    pub fn local_path(&self) -> Option<&std::path::Path> {
        match &self.storage {
            BamStorage::Local(path) => Some(path),
            BamStorage::Remote { .. } => None,
        }
    }

    pub fn remote_data(&self) -> Option<&[u8]> {
        match &self.storage {
            BamStorage::Local(_) => None,
            BamStorage::Remote { data, .. } => Some(data),
        }
    }

    pub fn remote_index_data(&self) -> Option<&[u8]> {
        match &self.storage {
            BamStorage::Local(_) => None,
            BamStorage::Remote { index_data, .. } => index_data.as_deref(),
        }
    }
}

/// Return candidate index URIs for a BAM location.
pub fn index_uri_candidates(bam_uri: &str) -> Vec<String> {
    if looks_like_uri(bam_uri) {
        vec![
            format!("{bam_uri}.bai"),
            replace_suffix(bam_uri, ".bam", ".bai"),
        ]
    } else {
        let path = std::path::Path::new(bam_uri);
        let mut candidates = Vec::new();
        let mut with_bam_bai = path.to_path_buf();
        with_bam_bai.set_extension("bam.bai");
        candidates.push(with_bam_bai.display().to_string());

        let mut with_bai = path.to_path_buf();
        with_bai.set_extension("bai");
        candidates.push(with_bai.display().to_string());
        candidates
    }
}

/// Read a BAM (and optional `.bai`) from a local path or cloud URI.
pub fn open_bam(uri: &str) -> Result<BamSource, IoError> {
    if let Some(path) = local_path(uri) {
        return Ok(BamSource {
            uri: uri.to_string(),
            storage: BamStorage::Local(path),
        });
    }

    let data = read_bytes(uri)?;
    let (index_data, index_uri) = load_first_available_index(uri)?;

    Ok(BamSource {
        uri: uri.to_string(),
        storage: BamStorage::Remote {
            uri: uri.to_string(),
            data: std::sync::Arc::from(data),
            index_data: index_data.map(std::sync::Arc::from),
            index_uri,
        },
    })
}

/// Read an object from a local path or cloud URI.
pub fn read_bytes(uri: &str) -> Result<Vec<u8>, IoError> {
    if looks_like_uri(uri) {
        read_object_store(uri)
    } else {
        Ok(std::fs::read(uri)?)
    }
}

/// Return true when an index appears to exist for `bam_uri`.
pub fn has_index(source: &BamSource) -> bool {
    if let Some(index_data) = source.remote_index_data() {
        return !index_data.is_empty();
    }

    let bam_uri = &source.uri;
    index_uri_candidates(bam_uri).iter().any(|candidate| {
        if looks_like_uri(candidate) {
            object_exists(candidate).unwrap_or(false)
        } else {
            std::path::Path::new(candidate).exists()
        }
    })
}

fn load_first_available_index(bam_uri: &str) -> Result<(Option<Vec<u8>>, Option<String>), IoError> {
    for candidate in index_uri_candidates(bam_uri) {
        let bytes = if looks_like_uri(&candidate) {
            match read_object_store(&candidate) {
                Ok(bytes) => bytes,
                Err(IoError::ObjectStore(object_store::Error::NotFound { .. })) => continue,
                Err(IoError::ObjectStore(object_store::Error::Generic { store, source }))
                    if source.to_string().contains("404") =>
                {
                    let _ = store;
                    continue;
                }
                Err(err) => return Err(err),
            }
        } else if std::path::Path::new(&candidate).exists() {
            std::fs::read(&candidate)?
        } else {
            continue;
        };

        return Ok((Some(bytes), Some(candidate)));
    }

    Ok((None, None))
}

fn object_exists(uri: &str) -> Result<bool, IoError> {
    if !looks_like_uri(uri) {
        return Ok(std::path::Path::new(uri).exists());
    }

    runtime::block_on(async {
        let (store, path) = parse_store(uri)?;
        match store.head(&path).await {
            Ok(_) => Ok(true),
            Err(object_store::Error::NotFound { .. }) => Ok(false),
            Err(err) => Err(err.into()),
        }
    })
}

fn read_object_store(uri: &str) -> Result<Vec<u8>, IoError> {
    runtime::block_on(async {
        let (store, path) = parse_store(uri)?;
        let result = store.get(&path).await?;
        Ok(result.bytes().await?.to_vec())
    })
}

fn parse_store(
    uri: &str,
) -> Result<(Box<dyn object_store::ObjectStore>, object_store::path::Path), IoError> {
    let url = url::Url::parse(uri).map_err(|err| IoError::InvalidUri(err.to_string()))?;
    object_store::parse_url_opts(&url, std::iter::empty::<(&str, String)>())
        .map_err(IoError::from)
}

fn looks_like_uri(value: &str) -> bool {
    value.contains("://")
}

fn local_path(uri: &str) -> Option<std::path::PathBuf> {
    if looks_like_uri(uri) {
        if let Some(rest) = uri.strip_prefix("file://") {
            return Some(std::path::PathBuf::from(rest));
        }
        return None;
    }
    Some(std::path::PathBuf::from(uri))
}

fn replace_suffix(value: &str, from: &str, to: &str) -> String {
    if let Some(prefix) = value.strip_suffix(from) {
        format!("{prefix}{to}")
    } else {
        format!("{value}{to}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_index_uri_candidates_for_cloud_paths() {
        let candidates = index_uri_candidates("s3://bucket/path/sample.bam");
        assert_eq!(
            candidates,
            vec![
                "s3://bucket/path/sample.bam.bai".to_string(),
                "s3://bucket/path/sample.bai".to_string(),
            ]
        );
    }

    #[test]
    fn reads_local_fixture_via_file_uri() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/data/tiny.bam");
        if !path.exists() {
            return;
        }

        let uri = format!("file://{}", path.display());
        let source = open_bam(&uri).expect("open file:// bam");
        assert!(source.local_path().is_some());
        assert_eq!(source.uri, uri);
    }
}