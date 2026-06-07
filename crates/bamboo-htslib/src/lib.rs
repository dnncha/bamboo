//! htslib-backed operations for Bamboo.
//!
//! BAM/VCF read paths use `bamboo-noodles`; pileup and future write paths use htslib.

#[cfg(feature = "htslib")]
mod error;
#[cfg(feature = "htslib")]
mod pileup;

#[cfg(feature = "htslib")]
pub use error::HtslibError;
#[cfg(feature = "htslib")]
pub use pileup::{pileup_region, PileupColumn, PileupRead};

/// Placeholder marker for future htslib integration.
pub const PHASE: &str = "phase-2-cram-via-noodles";

/// Returns whether htslib-backed features are compiled in.
pub fn is_available() -> bool {
    cfg!(feature = "htslib")
}

/// Returns whether pileup is available through htslib.
pub fn pileup_available() -> bool {
    cfg!(feature = "htslib")
}

/// Recommended backend for alignment I/O today.
pub fn primary_backend() -> &'static str {
    "noodles"
}