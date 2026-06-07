use crate::error::NoodlesError;
use bamboo_core::FetchRegion;
use noodles::fasta as fasta;
use noodles::fasta::record::{Definition, Sequence as FastaSequence};
use noodles::fasta::repository::adapters::IndexedReader;
use noodles::sam::Header;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

/// Build or reuse a CRAM reference repository tuned for the scan shape.
pub fn reference_repository_for_scan(
    header: &Header,
    reference_fasta: Option<&str>,
    region: Option<&FetchRegion>,
) -> Result<fasta::Repository, NoodlesError> {
    if let Some(path) = reference_fasta {
        return cached_fasta_repository(path, region);
    }
    Ok(reference_repository_from_header(header))
}

fn cached_fasta_repository(
    path: &str,
    region: Option<&FetchRegion>,
) -> Result<fasta::Repository, NoodlesError> {
    static CACHE: OnceLock<Mutex<HashMap<String, fasta::Repository>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));

    let cache_key = if let Some(region) = region {
        format!("{path}#{}", region.reference_name)
    } else {
        path.to_string()
    };

    let mut guard = cache.lock().expect("FASTA repository cache lock");
    if let Some(repository) = guard.get(&cache_key) {
        return Ok(repository.clone());
    }

    let repository = if fasta_index_path(path).exists() {
        reference_repository_from_indexed_fasta(path)?
    } else if let Some(region) = region {
        reference_repository_from_fasta_contig(path, &region.reference_name)?
    } else {
        reference_repository_from_fasta(path)?
    };

    guard.insert(cache_key, repository.clone());
    Ok(repository)
}

pub fn fasta_index_path(path: &str) -> PathBuf {
    let mut candidate = std::ffi::OsString::from(path);
    candidate.push(".fai");
    PathBuf::from(candidate)
}

fn reference_repository_from_indexed_fasta(path: &str) -> Result<fasta::Repository, NoodlesError> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
    Ok(fasta::Repository::new(IndexedReader::new(reader)))
}

fn reference_repository_from_fasta(path: &str) -> Result<fasta::Repository, NoodlesError> {
    let mut reader = fasta::io::reader::Builder::default()
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
    let records = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .map_err(NoodlesError::from)?;
    Ok(fasta::Repository::new(records))
}

fn reference_repository_from_fasta_contig(
    path: &str,
    contig: &str,
) -> Result<fasta::Repository, NoodlesError> {
    let mut reader = fasta::io::reader::Builder::default()
        .build_from_path(path)
        .map_err(NoodlesError::from)?;

    for result in reader.records() {
        let record = result.map_err(NoodlesError::from)?;
        if record.name() == contig.as_bytes() {
            return Ok(fasta::Repository::new(vec![record]));
        }
    }

    Err(NoodlesError::Message(format!(
        "reference contig '{contig}' not found in {path}"
    )))
}

fn reference_repository_from_header(header: &Header) -> fasta::Repository {
    let records = header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| {
            fasta::Record::new(
                Definition::new(name.to_string(), None),
                FastaSequence::from(vec![b'N'; reference.length().get()]),
            )
        })
        .collect::<Vec<_>>();
    fasta::Repository::new(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{tiny_fasta_path, write_tiny_fasta, write_tiny_fasta_index};
    use crate::header_util::header_from_references;
    use bamboo_core::FetchRegion;
    use tempfile::tempdir;

    #[test]
    fn header_repository_materializes_all_reference_sequences() {
        let header = header_from_references(&[
            ("chr1".to_string(), 1_000),
            ("chr2".to_string(), 2_000),
        ])
        .expect("header");

        let repository = reference_repository_from_header(&header);
        assert_eq!(repository.len(), 0);

        let chr1 = repository
            .get(b"chr1")
            .transpose()
            .expect("lookup result")
            .expect("chr1 sequence");
        assert_eq!(chr1.len(), 1_000);

        let chr2 = repository
            .get(b"chr2")
            .transpose()
            .expect("lookup result")
            .expect("chr2 sequence");
        assert_eq!(chr2.len(), 2_000);
    }

    #[test]
    fn indexed_fasta_repository_resolves_contig_without_loading_entire_file() {
        let dir = tempdir().unwrap();
        let fasta_path = tiny_fasta_path(dir.path());
        write_tiny_fasta(&fasta_path).unwrap();
        write_tiny_fasta_index(&fasta_path).unwrap();

        let repository =
            reference_repository_from_indexed_fasta(fasta_path.to_str().unwrap()).unwrap();
        assert_eq!(repository.len(), 0);

        let sequence = repository
            .get(b"chr1")
            .transpose()
            .expect("lookup result")
            .expect("chr1 sequence");
        assert_eq!(sequence.len(), 1_000);
    }
}