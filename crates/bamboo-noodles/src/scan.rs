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
    let mut table = BamTable::new(options.columns.clone(), options.tags.clone());
    for record in reader.open_stream(options.clone())? {
        record?.append_to_table(&mut table, &options);
    }
    Ok(table)
}

#[allow(dead_code)]
pub fn records_to_table(records: Vec<AlignedRecord>, options: &BamScanOptions) -> BamTable {
    let mut table = BamTable::new(options.columns.clone(), options.tags.clone());
    for record in records {
        record.append_to_table(&mut table, options);
    }
    table
}