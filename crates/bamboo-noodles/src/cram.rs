use crate::error::NoodlesError;
use crate::record::AlignedRecord;
use bamboo_core::BamScanOptions;
use noodles::cram as cram;
use noodles::fasta as fasta;
use noodles::fasta::record::{Definition, Sequence as FastaSequence};
use noodles::sam::Header;
use std::path::Path;
use std::vec;

/// Record stream over a CRAM file.
///
/// Records are decoded up front via noodles' public `records` iterator. This keeps
/// the implementation simple while CRAM indexed query support matures.
pub struct CramRecordStream {
    records: vec::IntoIter<AlignedRecord>,
}

impl CramRecordStream {
    pub fn open(path: &str, options: BamScanOptions) -> Result<Self, NoodlesError> {
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(reference_repository_from_header(
                &read_header_from_path(path)?,
            ))
            .build_from_path(path)
            .map_err(NoodlesError::from)?;
        let header = reader.read_header().map_err(NoodlesError::from)?;

        let mut records = Vec::new();
        for result in reader.records(&header) {
            let cram_record = result.map_err(NoodlesError::from)?;
            let record_buf = cram_record
                .try_into_alignment_record(&header)
                .map_err(NoodlesError::from)?;
            let aligned = AlignedRecord::from_record_buf(&header, &record_buf, &options);
            if aligned.passes_filters(&options) {
                records.push(aligned);
            }
        }

        Ok(Self {
            records: records.into_iter(),
        })
    }
}

impl Iterator for CramRecordStream {
    type Item = Result<AlignedRecord, NoodlesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next().map(Ok)
    }
}

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

    pub fn open_stream(&self, options: BamScanOptions) -> Result<CramRecordStream, NoodlesError> {
        CramRecordStream::open(&self.path, options)
    }

    pub fn iter_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        self.open_stream(options.clone())?
            .collect::<Result<Vec<_>, _>>()
    }

    pub fn fetch_records(&self, options: &BamScanOptions) -> Result<Vec<AlignedRecord>, NoodlesError> {
        self.iter_records(options)
    }

    pub fn count_records(&self) -> Result<usize, NoodlesError> {
        let options = BamScanOptions::iteration_defaults();
        Ok(self.iter_records(&options)?.len())
    }
}

fn read_header_from_path(path: &str) -> Result<Header, NoodlesError> {
    let mut reader = cram::io::reader::Builder::default()
        .build_from_path(Path::new(path))
        .map_err(NoodlesError::from)?;
    reader.read_header().map_err(NoodlesError::from)
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
    use crate::fixtures::{tiny_bam_path, tiny_cram_path, write_tiny_bam, write_tiny_cram};
    use bamboo_core::{BamColumn, FetchRegion};
    use tempfile::tempdir;

    fn iteration_options() -> BamScanOptions {
        BamScanOptions::iteration_defaults()
    }

    #[test]
    fn reads_tiny_cram() {
        let dir = tempdir().unwrap();
        let path = tiny_cram_path(dir.path());
        write_tiny_cram(&path).unwrap();

        let reader = CramReader::open(path.to_str().unwrap()).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
        assert_eq!(reader.reference_names(), vec!["chr1", "chr2"]);
    }

    #[test]
    fn iterates_tiny_cram_records() {
        let dir = tempdir().unwrap();
        let cram_path = tiny_cram_path(dir.path());
        let bam_path = tiny_bam_path(dir.path());
        write_tiny_cram(&cram_path).unwrap();
        write_tiny_bam(&bam_path).unwrap();

        let cram_reader = CramReader::open(cram_path.to_str().unwrap()).unwrap();
        let cram_records = cram_reader.iter_records(&iteration_options()).unwrap();

        let bam_reader = crate::BamReader::open(bam_path.to_str().unwrap()).unwrap();
        let bam_records = bam_reader.iter_records(&iteration_options()).unwrap();

        assert_eq!(cram_records.len(), bam_records.len());
        for (cram, bam) in cram_records.iter().zip(bam_records.iter()) {
            assert_eq!(cram.query_name, bam.query_name);
            assert_eq!(cram.reference_name, bam.reference_name);
            assert_eq!(cram.reference_start, bam.reference_start);
            assert_eq!(cram.mapping_quality, bam.mapping_quality);
            assert_eq!(cram.cigar, bam.cigar);
        }
    }

    #[test]
    fn filters_cram_by_region() {
        let dir = tempdir().unwrap();
        let path = tiny_cram_path(dir.path());
        write_tiny_cram(&path).unwrap();

        let reader = CramReader::open(path.to_str().unwrap()).unwrap();
        let options = BamScanOptions {
            region: Some(
                FetchRegion::from_samtools_region("chr1:100-101").expect("valid region"),
            ),
            columns: vec![BamColumn::QueryName, BamColumn::ReferenceName, BamColumn::Position],
            ..Default::default()
        };
        let records = reader.iter_records(&options).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name.as_deref(), Some("read1"));
    }
}