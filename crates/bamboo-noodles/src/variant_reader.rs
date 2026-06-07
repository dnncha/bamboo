use crate::bcf::BcfReader;
use crate::error::NoodlesError;
use crate::vcf::VcfReader;
use bamboo_core::{VcfScanOptions, VcfTable};

/// Unified variant reader for VCF and BCF inputs.
pub enum VariantReader {
    Vcf(VcfReader),
    Bcf(BcfReader),
}

impl VariantReader {
    pub fn open(path: &str) -> Result<Self, NoodlesError> {
        if path
            .rsplit_once('.')
            .is_some_and(|(_, extension)| extension.eq_ignore_ascii_case("bcf"))
        {
            Ok(Self::Bcf(BcfReader::open(path)?))
        } else {
            Ok(Self::Vcf(VcfReader::open(path)?))
        }
    }

    pub fn path(&self) -> &str {
        match self {
            Self::Vcf(reader) => reader.path(),
            Self::Bcf(reader) => reader.path(),
        }
    }

    pub fn reference_names(&self) -> Vec<String> {
        match self {
            Self::Vcf(reader) => reader.reference_names(),
            Self::Bcf(reader) => reader.reference_names(),
        }
    }

    pub fn has_index(&self) -> bool {
        match self {
            Self::Vcf(reader) => reader.has_index(),
            Self::Bcf(reader) => reader.has_index(),
        }
    }

    pub fn count_records(&self) -> Result<usize, NoodlesError> {
        match self {
            Self::Vcf(reader) => reader.count_records(),
            Self::Bcf(reader) => reader.count_records(),
        }
    }

    pub fn scan(&self, options: VcfScanOptions) -> Result<VcfTable, NoodlesError> {
        match self {
            Self::Vcf(reader) => reader.scan(options),
            Self::Bcf(reader) => reader.scan(options),
        }
    }
}