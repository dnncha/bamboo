use crate::error::NoodlesError;
use crate::record::AlignedRecord;
use bamboo_core::{BamScanOptions, FetchRegion};
use noodles::bam as bam;
use noodles::bgzf as bgzf;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use std::path::{Path, PathBuf};

/// High-level BAM reader backed by noodles.
pub struct BamReader {
    path: PathBuf,
    header: Header,
    data: Vec<u8>,
}

impl BamReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, NoodlesError> {
        let path = path.as_ref().to_path_buf();
        let data = std::fs::read(&path)?;
        let header = read_header(&data)?;
        Ok(Self { path, header, data })
    }

    pub fn path(&self) -> &Path {
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

    pub fn has_index(&self) -> bool {
        index_path(&self.path).exists()
    }

    pub fn count_records(&self) -> Result<usize, NoodlesError> {
        let mut reader = open_plain_reader(&self.data)?;
        let _ = reader.read_header().map_err(NoodlesError::from)?;
        Ok(reader.record_bufs(&self.header).count())
    }

    pub fn iter_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        if options.region.is_some() && index_path(&self.path).exists() {
            self.fetch_records(options)
        } else {
            self.scan_records(options)
        }
    }

    pub fn scan_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        let mut reader = open_plain_reader(&self.data)?;
        let _ = reader.read_header().map_err(NoodlesError::from)?;
        let mut records = Vec::new();

        for result in reader.record_bufs(&self.header) {
            let record = result?;
            let aligned = AlignedRecord::from_record_buf(&self.header, record, &options.tags);
            if aligned.passes_filters(options) {
                records.push(aligned);
            }
        }

        Ok(records)
    }

    pub fn fetch_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        let region = options
            .region
            .as_ref()
            .ok_or_else(|| NoodlesError::Message("fetch requires a region".to_string()))?;

        if !index_path(&self.path).exists() {
            return Err(NoodlesError::MissingIndex {
                path: index_path(&self.path).display().to_string(),
            });
        }

        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&self.path)
            .map_err(NoodlesError::from)?;

        let header = reader.read_header().map_err(NoodlesError::from)?;
        let samtools_region = region.to_samtools_region();
        let parsed_region: noodles::core::Region = samtools_region
            .parse()
            .map_err(|err: noodles::core::region::ParseError| {
                NoodlesError::Message(err.to_string())
            })?;

        let query = reader
            .query(&header, &parsed_region)
            .map_err(NoodlesError::from)?;

        let mut records = Vec::new();
        for result in query {
            let record = result.map_err(NoodlesError::from)?;
            let record_buf =
                RecordBuf::try_from_alignment_record(&header, &record).map_err(NoodlesError::from)?;
            let aligned = AlignedRecord::from_record_buf(&header, record_buf, &options.tags);
            if aligned.passes_filters(options) {
                records.push(aligned);
            }
        }

        Ok(records)
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

fn read_header(data: &[u8]) -> Result<Header, NoodlesError> {
    let mut reader = open_plain_reader(data)?;
    reader.read_header().map_err(NoodlesError::from)
}

fn open_plain_reader(data: &[u8]) -> Result<bam::io::Reader<bgzf::Reader<&[u8]>>, NoodlesError> {
    Ok(bam::io::Reader::new(data))
}

fn index_path(path: &Path) -> PathBuf {
    let mut index_path = path.to_path_buf();
    index_path.set_extension("bam.bai");
    if index_path.exists() {
        return index_path;
    }

    let mut alt = path.to_path_buf();
    alt.set_extension("bai");
    alt
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

        let reader = BamReader::open(&path).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
    }

    #[test]
    fn reads_tiny_bam_records() {
        let dir = tempdir().unwrap();
        let path = tiny_bam_path(dir.path());
        write_tiny_bam(&path).unwrap();

        let reader = BamReader::open(&path).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
        assert_eq!(reader.reference_names(), vec!["chr1", "chr2"]);
    }

    #[test]
    fn filters_by_mapq() {
        let dir = tempdir().unwrap();
        let path = tiny_bam_path(dir.path());
        write_tiny_bam(&path).unwrap();

        let reader = BamReader::open(&path).unwrap();
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

        let reader = BamReader::open(&path).unwrap();
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

        let reader = BamReader::open(&path).unwrap();
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
}