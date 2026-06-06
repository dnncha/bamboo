use crate::region::FetchRegion;
use crate::schema::DEFAULT_BAM_COLUMNS;

/// Fixed alignment columns supported by Bamboo scanners.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BamColumn {
    QueryName,
    Flag,
    ReferenceName,
    Position,
    MappingQuality,
    Cigar,
    MateReferenceName,
    MatePosition,
    TemplateLength,
    Sequence,
    Quality,
}

impl BamColumn {
    pub fn parse_name(name: &str) -> Option<Self> {
        match name {
            "qname" => Some(Self::QueryName),
            "flag" => Some(Self::Flag),
            "rname" => Some(Self::ReferenceName),
            "pos" => Some(Self::Position),
            "mapq" => Some(Self::MappingQuality),
            "cigar" => Some(Self::Cigar),
            "rnext" => Some(Self::MateReferenceName),
            "pnext" => Some(Self::MatePosition),
            "tlen" => Some(Self::TemplateLength),
            "seq" => Some(Self::Sequence),
            "qual" => Some(Self::Quality),
            _ => None,
        }
    }
}

/// Query options for BAM scanning with projection and predicate pushdown.
#[derive(Debug, Clone)]
pub struct BamScanOptions {
    pub columns: Vec<BamColumn>,
    pub tags: Vec<String>,
    pub region: Option<FetchRegion>,
    pub min_mapq: Option<u8>,
    pub reference_name: Option<String>,
}

impl Default for BamScanOptions {
    fn default() -> Self {
        Self {
            columns: DEFAULT_BAM_COLUMNS.to_vec(),
            tags: Vec::new(),
            region: None,
            min_mapq: None,
            reference_name: None,
        }
    }
}

impl BamScanOptions {
    pub fn wants_column(&self, column: BamColumn) -> bool {
        self.columns.iter().any(|c| *c == column)
    }

    pub fn wants_tag(&self, tag: &str) -> bool {
        self.tags.iter().any(|t| t == tag)
    }
}