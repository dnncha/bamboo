use crate::error::NoodlesError;
use noodles::cram as cram;
use noodles::sam::Header;
use std::path::Path;

/// High-level CRAM reader backed by noodles.
pub struct CramReader {
    path: String,
    header: Header,
}

impl CramReader {
    pub fn open(path: &str) -> Result<Self, NoodlesError> {
        let header = read_header_from_path(path)?;
        Ok(Self {
            path: path.to_string(),
            header,
        })
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn reference_names(&self) -> Vec<String> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect()
    }

    pub fn reference_lengths(&self) -> Vec<u32> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(_, reference)| reference.length().get() as u32)
            .collect()
    }

    pub fn count_records(&self) -> Result<usize, NoodlesError> {
        let mut reader = cram::io::reader::Builder::default()
            .build_from_path(&self.path)
            .map_err(NoodlesError::from)?;
        reader.read_header().map_err(NoodlesError::from)?;

        let mut total = 0usize;
        while let Some(container) = reader.read_data_container().map_err(NoodlesError::from)? {
            for slice in container.slices() {
                let records = slice
                    .records(container.compression_header())
                    .map_err(NoodlesError::from)?;
                total += records.len();
            }
        }

        Ok(total)
    }
}

fn read_header_from_path(path: &str) -> Result<Header, NoodlesError> {
    let mut reader = cram::io::reader::Builder::default()
        .build_from_path(Path::new(path))
        .map_err(NoodlesError::from)?;
    reader.read_header().map_err(NoodlesError::from)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{tiny_cram_path, write_tiny_cram};
    use tempfile::tempdir;

    #[test]
    fn reads_tiny_cram() {
        let dir = tempdir().unwrap();
        let path = tiny_cram_path(dir.path());
        write_tiny_cram(&path).unwrap();

        let reader = CramReader::open(path.to_str().unwrap()).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
        assert_eq!(reader.reference_names(), vec!["chr1", "chr2"]);
    }
}