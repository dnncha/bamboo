//! Core types shared across Bamboo backends.

mod options;
mod region;
mod schema;
mod table;
mod vcf_table;

pub use options::{BamColumn, BamScanOptions, ITERATION_BAM_COLUMNS};
pub use region::{FetchRegion, RegionParseError};
pub use schema::{DEFAULT_BAM_COLUMNS, bam_column_name};
pub use table::{BamTable, TagColumn, TagValue};
pub use vcf_table::{DEFAULT_VCF_COLUMNS, VcfColumn, VcfScanOptions, VcfTable};