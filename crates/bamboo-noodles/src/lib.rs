//! Noodles-backed BAM reader and scanner.

pub mod fixtures;
mod error;
mod reader;
mod record;
mod scan;

pub use error::NoodlesError;
pub use reader::BamReader;
pub use record::AlignedRecord;
pub use scan::scan_bam;