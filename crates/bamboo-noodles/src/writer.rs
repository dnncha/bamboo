use crate::error::NoodlesError;
use crate::header_util::header_from_references;
use crate::record::AlignedRecord;
use noodles::bam as bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write;
use std::path::Path;

/// High-level BAM writer backed by noodles.
pub struct BamWriter {
    writer: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
    header: Header,
    path: String,
    finished: bool,
}

impl BamWriter {
    /// Create a new BAM at a local filesystem path.
    pub fn create(path: &str, header: Header) -> Result<Self, NoodlesError> {
        if path.contains("://") {
            return Err(NoodlesError::Message(
                "BAM writing only supports local paths (cloud write is not implemented)".to_string(),
            ));
        }

        let mut writer = bam::io::writer::Builder::default()
            .build_from_path(path)
            .map_err(NoodlesError::from)?;
        writer
            .write_header(&header)
            .map_err(NoodlesError::from)?;

        Ok(Self {
            writer,
            header,
            path: path.to_string(),
            finished: false,
        })
    }

    /// Convenience constructor from reference dictionaries.
    pub fn create_with_references(path: &str, refs: &[(String, u32)]) -> Result<Self, NoodlesError> {
        Self::create(path, header_from_references(refs)?)
    }

    /// Create a writer using the header from an open reader (pysam-style template).
    pub fn create_from_reader(path: &str, reader: &crate::reader::BamReader) -> Result<Self, NoodlesError> {
        Self::create(path, reader.header().clone())
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn write_record(&mut self, record: &AlignedRecord) -> Result<(), NoodlesError> {
        let record_buf = record.to_record_buf(&self.header)?;
        self.writer
            .write_alignment_record(&self.header, &record_buf)
            .map_err(NoodlesError::from)
    }

    pub fn finish(&mut self) -> Result<(), NoodlesError> {
        if !self.finished {
            self.writer.try_finish().map_err(NoodlesError::from)?;
            self.finished = true;
        }
        Ok(())
    }

    /// Build a `.bai` index for a finished BAM on disk.
    pub fn write_index(path: &str) -> Result<(), NoodlesError> {
        if path.contains("://") {
            return Err(NoodlesError::Message(
                "BAM index writing only supports local paths".to_string(),
            ));
        }

        let bam_path = Path::new(path);
        let index = bam::fs::index(bam_path).map_err(NoodlesError::from)?;
        let index_path = bam_path.with_extension("bam.bai");
        bam::bai::fs::write(index_path, &index).map_err(NoodlesError::from)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{tiny_bam_path, write_tiny_bam};
    use crate::reader::BamReader;
    use bamboo_core::{BamColumn, BamScanOptions};
    use tempfile::tempdir;

    #[test]
    fn round_trips_tiny_bam_through_writer() {
        let dir = tempdir().unwrap();
        let source_path = tiny_bam_path(dir.path());
        write_tiny_bam(&source_path).unwrap();

        let reader = BamReader::open(source_path.to_str().unwrap()).unwrap();
        let refs: Vec<_> = reader
            .reference_names()
            .into_iter()
            .zip(reader.reference_lengths())
            .collect();

        let out_path = dir.path().join("copy.bam");
        let mut writer =
            BamWriter::create_with_references(out_path.to_str().unwrap(), &refs).unwrap();

        let options = BamScanOptions {
            columns: vec![
                BamColumn::QueryName,
                BamColumn::ReferenceName,
                BamColumn::Position,
                BamColumn::MappingQuality,
                BamColumn::Cigar,
                BamColumn::Sequence,
                BamColumn::Quality,
            ],
            ..Default::default()
        };
        for record in reader.scan_records(&options).unwrap() {
            writer.write_record(&record).unwrap();
        }
        writer.finish().unwrap();

        let copied = BamReader::open(out_path.to_str().unwrap()).unwrap();
        assert_eq!(copied.count_records().unwrap(), 2);
        assert_eq!(copied.reference_names(), vec!["chr1", "chr2"]);
    }
}