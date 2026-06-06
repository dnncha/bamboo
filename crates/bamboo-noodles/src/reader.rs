use crate::error::NoodlesError;
use crate::record::AlignedRecord;
use crate::stream::{BamRecordStream, count_records};
use bamboo_core::{BamScanOptions, FetchRegion};
use bamboo_io::BamSource;
use noodles::bam as bam;
use noodles::bgzf as bgzf;
use noodles::sam::Header;
use std::io::Cursor;

/// High-level BAM reader backed by noodles.
pub struct BamReader {
    source: BamSource,
    header: Header,
}

impl BamReader {
    /// Open a BAM from a local path or cloud URI (`s3://`, `gs://`, `https://`, `file://`).
    pub fn open(uri: &str) -> Result<Self, NoodlesError> {
        let source = bamboo_io::open_bam(uri).map_err(NoodlesError::from)?;
        let header = read_header(&source)?;
        Ok(Self { source, header })
    }

    pub fn uri(&self) -> &str {
        &self.source.uri
    }

    pub fn source(&self) -> &BamSource {
        &self.source
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

    pub fn has_index(&self) -> bool {
        bamboo_io::has_index(&self.source)
    }

    pub fn count_records(&self) -> Result<usize, NoodlesError> {
        count_records(&self.source, &self.header)
    }

    pub fn open_stream(&self, options: BamScanOptions) -> Result<BamRecordStream, NoodlesError> {
        BamRecordStream::open(&self.source, &self.header, options)
    }

    pub fn iter_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        self.open_stream(options.clone())?
            .collect::<Result<Vec<_>, _>>()
    }

    pub fn scan_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        self.iter_records(options)
    }

    pub fn fetch_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        if !self.has_index() {
            return Err(NoodlesError::MissingIndex {
                path: bamboo_io::index_uri_candidates(&self.source.uri).join(", "),
            });
        }
        self.iter_records(options)
    }

    pub fn fetch_region(
        &self,
        region: &FetchRegion,
        tags: &[String],
    ) -> Result<Vec<AlignedRecord>, NoodlesError> {
        let options = BamScanOptions {
            region: Some(region.clone()),
            tags: tags.to_vec(),
            ..BamScanOptions::default()
        };
        self.fetch_records(&options)
    }
}

fn read_header(source: &BamSource) -> Result<Header, NoodlesError> {
    match &source.storage {
        bamboo_io::BamStorage::Local(path) => {
            let mut reader = bam::io::reader::Builder::default()
                .build_from_path(path)
                .map_err(NoodlesError::from)?;
            reader.read_header().map_err(NoodlesError::from)
        }
        bamboo_io::BamStorage::Remote { data, .. } => {
            let mut reader = open_plain_reader(data)?;
            reader.read_header().map_err(NoodlesError::from)
        }
    }
}

fn open_plain_reader(
    data: &std::sync::Arc<[u8]>,
) -> Result<bam::io::Reader<bgzf::Reader<Cursor<std::sync::Arc<[u8]>>>>, NoodlesError> {
    Ok(bam::io::Reader::new(Cursor::new(data.clone())))
}

#[cfg(test)]
mod tests {
    use crate::fixtures::{tiny_bam_path, write_tiny_bam, write_tiny_bam_index};
    use super::BamReader;
    use bamboo_core::{BamColumn, BamScanOptions, FetchRegion};
    use noodles::bam as bam;
    use tempfile::tempdir;

    #[test]
    fn fixture_roundtrip_writes_two_records() {
        let dir = tempdir().unwrap();
        let path = tiny_bam_path(dir.path());
        write_tiny_bam(&path).unwrap();

        let data = std::fs::read(&path).unwrap();
        let mut reader = bam::io::Reader::new(data.as_slice());
        let header = reader.read_header().unwrap();
        let records: Vec<_> = reader
            .record_bufs(&header)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(records.len(), 2);

        let reader = BamReader::open(path.to_str().unwrap()).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
    }

    #[test]
    fn reads_tiny_bam_records() {
        let dir = tempdir().unwrap();
        let path = tiny_bam_path(dir.path());
        write_tiny_bam(&path).unwrap();

        let reader = BamReader::open(path.to_str().unwrap()).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
        assert_eq!(reader.reference_names(), vec!["chr1", "chr2"]);
    }

    #[test]
    fn filters_by_mapq() {
        let dir = tempdir().unwrap();
        let path = tiny_bam_path(dir.path());
        write_tiny_bam(&path).unwrap();

        let reader = BamReader::open(path.to_str().unwrap()).unwrap();
        let options = BamScanOptions {
            min_mapq: Some(30),
            columns: vec![BamColumn::QueryName, BamColumn::MappingQuality],
            ..Default::default()
        };
        let records = reader.scan_records(&options).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name.as_deref(), Some("read1"));
    }

    #[test]
    fn filters_by_region_without_index() {
        let dir = tempdir().unwrap();
        let path = tiny_bam_path(dir.path());
        write_tiny_bam(&path).unwrap();

        let reader = BamReader::open(path.to_str().unwrap()).unwrap();
        let options = BamScanOptions {
            region: Some(FetchRegion {
                reference_name: "chr1".to_string(),
                start: Some(50),
                end: Some(150),
            }),
            columns: vec![BamColumn::QueryName],
            ..Default::default()
        };
        let records = reader.scan_records(&options).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name.as_deref(), Some("read1"));
    }

    #[test]
    fn fetches_by_index_when_available() {
        let dir = tempdir().unwrap();
        let path = tiny_bam_path(dir.path());
        write_tiny_bam(&path).unwrap();
        write_tiny_bam_index(&path).unwrap();

        let reader = BamReader::open(path.to_str().unwrap()).unwrap();
        assert!(reader.has_index());

        let options = BamScanOptions {
            region: Some(FetchRegion {
                reference_name: "chr1".to_string(),
                start: Some(50),
                end: Some(150),
            }),
            columns: vec![BamColumn::QueryName],
            ..Default::default()
        };
        let records = reader.fetch_records(&options).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name.as_deref(), Some("read1"));
    }

    #[test]
    fn opens_file_uri() {
        let dir = tempdir().unwrap();
        let path = tiny_bam_path(dir.path());
        write_tiny_bam(&path).unwrap();
        write_tiny_bam_index(&path).unwrap();

        let uri = format!("file://{}", path.display());
        let reader = BamReader::open(&uri).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
        assert!(reader.has_index());
    }
}