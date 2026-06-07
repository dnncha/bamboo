use crate::error::NoodlesError;
use crate::record::AlignedRecord;
use bamboo_core::BamScanOptions;
use noodles::cram as cram;
use noodles::fasta as fasta;
use noodles::fasta::record::{Definition, Sequence as FastaSequence};
use noodles::sam::Header;
use std::path::{Path, PathBuf};
use std::vec;

/// Record stream over a CRAM file.
pub struct CramRecordStream {
    records: vec::IntoIter<AlignedRecord>,
}

impl CramRecordStream {
    pub fn open(
        path: &str,
        options: BamScanOptions,
        reference_fasta: Option<&str>,
    ) -> Result<Self, NoodlesError> {
        let header = read_header_from_path(path)?;
        let repository = reference_repository(&header, reference_fasta)?;

        let mut records = Vec::new();
        if let Some(region) = &options.region {
            let index_path = cram_index_path(path);
            if index_path.exists() {
                collect_indexed_region_records(path, &header, &repository, region, &options, &mut records)?;
                return Ok(Self {
                    records: records.into_iter(),
                });
            }
        }

        collect_linear_records(path, &header, &repository, &options, &mut records)?;
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
    reference_fasta: Option<String>,
}

impl CramReader {
    pub fn open(path: &str) -> Result<Self, NoodlesError> {
        Self::open_with_reference(path, None)
    }

    pub fn open_with_reference(
        path: &str,
        reference_fasta: Option<&str>,
    ) -> Result<Self, NoodlesError> {
        let header = read_header_from_path(path)?;
        Ok(Self {
            path: path.to_string(),
            header,
            reference_fasta: reference_fasta.map(str::to_string),
        })
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn reference_fasta(&self) -> Option<&str> {
        self.reference_fasta.as_deref()
    }

    pub fn has_index(&self) -> bool {
        cram_index_path(&self.path).exists()
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
        CramRecordStream::open(&self.path, options, self.reference_fasta.as_deref())
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

    pub fn scan(&self, options: BamScanOptions) -> Result<bamboo_core::BamTable, NoodlesError> {
        crate::scan::scan_cram_reader(self, options)
    }
}

pub(crate) fn cram_index_path(path: &str) -> PathBuf {
    let mut name = PathBuf::from(path).into_os_string();
    name.push(".crai");
    PathBuf::from(name)
}

fn collect_indexed_region_records(
    path: &str,
    _header: &Header,
    repository: &fasta::Repository,
    region: &bamboo_core::FetchRegion,
    options: &BamScanOptions,
    records: &mut Vec<AlignedRecord>,
) -> Result<(), NoodlesError> {
    let parsed_region: noodles::core::Region = region
        .to_samtools_region()
        .parse()
        .map_err(|err: noodles::core::region::ParseError| NoodlesError::Message(err.to_string()))?;

    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
    let header = reader.read_header().map_err(NoodlesError::from)?;

    let query = reader
        .query(&header, &parsed_region)
        .map_err(NoodlesError::from)?;
    for result in query {
        push_aligned_record(records, result.map_err(NoodlesError::from)?, &header, options)?;
    }

    Ok(())
}

fn collect_linear_records(
    path: &str,
    _header: &Header,
    repository: &fasta::Repository,
    options: &BamScanOptions,
    records: &mut Vec<AlignedRecord>,
) -> Result<(), NoodlesError> {
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
    let header = reader.read_header().map_err(NoodlesError::from)?;

    for result in reader.records(&header) {
        push_aligned_record(records, result.map_err(NoodlesError::from)?, &header, options)?;
    }

    Ok(())
}

fn push_aligned_record(
    records: &mut Vec<AlignedRecord>,
    cram_record: cram::Record,
    header: &Header,
    options: &BamScanOptions,
) -> Result<(), NoodlesError> {
    let record_buf = cram_record
        .try_into_alignment_record(header)
        .map_err(NoodlesError::from)?;
    let aligned = AlignedRecord::from_record_buf(header, &record_buf, options);
    if aligned.passes_filters(options) {
        records.push(aligned);
    }
    Ok(())
}

pub(crate) fn read_header_from_path(path: &str) -> Result<Header, NoodlesError> {
    let mut reader = cram::io::reader::Builder::default()
        .build_from_path(Path::new(path))
        .map_err(NoodlesError::from)?;
    reader.read_header().map_err(NoodlesError::from)
}

pub(crate) fn reference_repository(
    header: &Header,
    reference_fasta: Option<&str>,
) -> Result<fasta::Repository, NoodlesError> {
    if let Some(path) = reference_fasta {
        return reference_repository_from_fasta(path);
    }
    Ok(reference_repository_from_header(header))
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
    use crate::fixtures::{
        tiny_bam_path, tiny_cram_index_path, tiny_cram_path, tiny_fasta_path, write_tiny_bam,
        write_tiny_cram, write_tiny_cram_index, write_tiny_fasta,
    };
    use bamboo_core::{BamColumn, FetchRegion};
    use tempfile::tempdir;

    fn iteration_options() -> BamScanOptions {
        BamScanOptions::iteration_defaults()
    }

    fn region_options(region: &str) -> BamScanOptions {
        BamScanOptions {
            region: Some(FetchRegion::from_samtools_region(region).expect("valid region")),
            columns: vec![BamColumn::QueryName, BamColumn::ReferenceName, BamColumn::Position],
            ..Default::default()
        }
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
        let records = reader.iter_records(&region_options("chr1:100-101")).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name.as_deref(), Some("read1"));
    }

    #[test]
    fn indexed_fetch_matches_linear_scan() {
        let dir = tempdir().unwrap();
        let path = tiny_cram_path(dir.path());
        write_tiny_cram(&path).unwrap();
        write_tiny_cram_index(&path).unwrap();

        let reader = CramReader::open(path.to_str().unwrap()).unwrap();
        assert!(reader.has_index());

        let options = region_options("chr1:100-101");
        let indexed = reader.iter_records(&options).unwrap();

        let index_path = tiny_cram_index_path(dir.path());
        assert!(index_path.exists());

        let linear = CramRecordStream::open(
            path.to_str().unwrap(),
            BamScanOptions {
                region: None,
                ..options.clone()
            },
            None,
        )
        .unwrap()
        .filter_map(Result::ok)
        .filter(|record| record.passes_filters(&options))
        .collect::<Vec<_>>();

        assert_eq!(indexed, linear);
    }

    #[test]
    fn scans_tiny_cram_to_columnar_table() {
        let dir = tempdir().unwrap();
        let cram_path = tiny_cram_path(dir.path());
        let bam_path = tiny_bam_path(dir.path());
        write_tiny_cram(&cram_path).unwrap();
        write_tiny_bam(&bam_path).unwrap();

        let options = BamScanOptions {
            columns: vec![BamColumn::QueryName, BamColumn::Position, BamColumn::MappingQuality],
            ..Default::default()
        };

        let cram_table = crate::scan::scan_cram(cram_path.to_str().unwrap(), options.clone(), None).unwrap();
        let bam_reader = crate::BamReader::open(bam_path.to_str().unwrap()).unwrap();
        let bam_table = crate::scan::scan_reader(&bam_reader, options).unwrap();

        assert_eq!(cram_table.len(), bam_table.len());
        assert_eq!(cram_table.qname, bam_table.qname);
        assert_eq!(cram_table.pos, bam_table.pos);
        assert_eq!(cram_table.mapq, bam_table.mapq);
    }

    #[test]
    fn reads_with_external_fasta_reference() {
        let dir = tempdir().unwrap();
        let cram_path = tiny_cram_path(dir.path());
        let fasta_path = tiny_fasta_path(dir.path());
        write_tiny_cram(&cram_path).unwrap();
        write_tiny_fasta(&fasta_path).unwrap();

        let reader =
            CramReader::open_with_reference(cram_path.to_str().unwrap(), Some(fasta_path.to_str().unwrap()))
                .unwrap();
        let records = reader.iter_records(&iteration_options()).unwrap();
        assert_eq!(records.len(), 2);
    }
}