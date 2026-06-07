//! Noodles-backed BAM reader and scanner.

pub mod fixtures;
mod bcf;
mod columnar;
mod cram;
mod lazy_fetch;
mod error;
mod header_util;
mod reader;
mod record;
mod scan;
mod stream;
mod variant_reader;
mod vcf;
mod writer;

pub use cram::{CramReader, CramRecordStream};
pub use error::NoodlesError;
pub use header_util::header_from_references;
pub use reader::BamReader;
pub use stream::BamRecordStream;
pub use record::AlignedRecord;
pub use scan::{scan_bam, scan_cram, scan_cram_reader, scan_reader};
pub use bcf::{scan_bcf, BcfReader};
pub use variant_reader::VariantReader;
pub use vcf::{scan_vcf, VcfReader};
pub use writer::BamWriter;