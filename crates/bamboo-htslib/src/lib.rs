//! htslib-backed operations for Bamboo.
//!
//! Phase 2 CRAM decoding with reference genomes and pileup remain here.
//! BAM/VCF read paths are implemented via `bamboo-noodles` (noodles).

/// Placeholder marker for future htslib integration.
pub const PHASE: &str = "phase-2-noodles-primary";

/// Returns whether htslib-backed features are compiled in.
pub fn is_available() -> bool {
    false
}

/// Recommended backend for alignment I/O today.
pub fn primary_backend() -> &'static str {
    "noodles"
}