//! Noodles-backed BAM reader and scanner.

pub mod fixtures;
mod error;
mod header_util;
mod reader;
mod record;
mod scan;
mod writer;

pub use error::NoodlesError;
pub use header_util::header_from_references;
pub use reader::BamReader;
pub use record::AlignedRecord;
pub use scan::scan_bam;
pub use writer::BamWriter;