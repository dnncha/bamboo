//! Noodles-backed BAM reader and scanner.

pub mod fixtures;
mod columnar;
mod lazy_fetch;
mod error;
mod header_util;
mod reader;
mod record;
mod scan;
mod stream;
mod vcf;
mod writer;

pub use error::NoodlesError;
pub use header_util::header_from_references;
pub use reader::BamReader;
pub use stream::BamRecordStream;
pub use record::AlignedRecord;
pub use scan::{scan_bam, scan_reader};
pub use vcf::{scan_vcf, VcfReader};
pub use writer::BamWriter;