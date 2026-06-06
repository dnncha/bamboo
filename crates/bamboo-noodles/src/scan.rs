use crate::columnar::scan_reader_columnar;
use crate::cram::CramReader;
use crate::error::NoodlesError;
use crate::reader::BamReader;
use crate::record::AlignedRecord;
use bamboo_core::{BamScanOptions, BamTable};

/// Scan a BAM file into a columnar `BamTable`.
pub fn scan_bam(path: &str, options: BamScanOptions) -> Result<BamTable, NoodlesError> {
    let reader = BamReader::open(path)?;
    scan_reader(&reader, options)
}

/// Scan an open reader into a columnar `BamTable` without reopening the source.
pub fn scan_reader(reader: &BamReader, options: BamScanOptions) -> Result<BamTable, NoodlesError> {
    scan_reader_columnar(reader.source(), reader.header(), options)
}

/// Scan a CRAM file into a columnar `BamTable`.
pub fn scan_cram(
    path: &str,
    options: BamScanOptions,
    reference_fasta: Option<&str>,
) -> Result<BamTable, NoodlesError> {
    let reader = CramReader::open_with_reference(path, reference_fasta)?;
    scan_cram_reader(&reader, options)
}

/// Scan an open CRAM reader into a columnar `BamTable`.
pub fn scan_cram_reader(reader: &CramReader, options: BamScanOptions) -> Result<BamTable, NoodlesError> {
    let records = reader.iter_records(&options)?;
    Ok(records_to_table(records, &options))
}

pub fn records_to_table(records: Vec<AlignedRecord>, options: &BamScanOptions) -> BamTable {
    let mut table = BamTable::new(options.columns.clone(), options.tags.clone());
    for record in records {
        record.append_to_table(&mut table, options);
    }
    table
}