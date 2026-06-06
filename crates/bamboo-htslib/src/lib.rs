//! htslib-backed operations for Bamboo.
//!
//! Phase 2 will add CRAM decoding, BAM/SAM/CRAM writing, and pileup here.
//! The MVP routes all BAM read paths through `bamboo-noodles`.

/// Placeholder marker for future htslib integration.
pub const PHASE: &str = "phase-2-stub";

/// Returns whether htslib-backed features are compiled in.
pub fn is_available() -> bool {
    false
}