use crate::error::NoodlesError;
use crate::record::AlignedRecord;
use bamboo_core::{BamScanOptions, FetchRegion};
use bamboo_io::BamSource;
use noodles::bam as bam;
use noodles::bgzf as bgzf;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
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
        let header = read_header(&source.data)?;
        Ok(Self { source, header })
    }

    pub fn uri(&self) -> &str {
        &self.source.uri
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
        bamboo_io::has_index(&self.source.uri, &self.source.index_data)
    }

    pub fn count_records(&self) -> Result<usize, NoodlesError> {
        let mut reader = open_plain_reader(&self.source.data)?;
        let _ = reader.read_header().map_err(NoodlesError::from)?;
        Ok(reader.record_bufs(&self.header).count())
    }

    pub fn iter_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        if options.region.is_some() && self.has_index() {
            self.fetch_records(options)
        } else {
            self.scan_records(options)
        }
    }

    pub fn scan_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        let mut reader = open_plain_reader(&self.source.data)?;
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

        let index_data = self
            .source
            .index_data
            .as_ref()
            .ok_or_else(|| NoodlesError::MissingIndex {
                path: self
                    .source
                    .index_uri
                    .clone()
                    .unwrap_or_else(|| bamboo_io::index_uri_candidates(&self.source.uri).join(", ")),
            })?;

        let mut index_reader = bam::bai::io::Reader::new(index_data.as_slice());
        let index = index_reader.read_index().map_err(NoodlesError::from)?;

        let mut reader = bam::io::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_reader(Cursor::new(self.source.data.as_slice()))
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