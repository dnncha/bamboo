use crate::columnar::scan_reader_columnar;
use crate::cram::CramReader;
use crate::error::NoodlesError;
use crate::reader::BamReader;

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
    crate::cram_columnar::scan_cram_columnar(path, options, reference_fasta)
}

/// Scan an open CRAM reader into a columnar `BamTable`.
pub fn scan_cram_reader(reader: &CramReader, options: BamScanOptions) -> Result<BamTable, NoodlesError> {
    crate::cram_columnar::scan_cram_columnar(
        reader.path(),
        options,
        reader.reference_fasta().as_deref(),
    )
}

