use crate::region::FetchRegion;

/// Fixed VCF columns supported by Bamboo scanners.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VcfColumn {
    Chrom,
    Pos,
    Id,
    Ref,
    Alt,
    Qual,
    Filter,
}

impl VcfColumn {
    pub fn parse_name(name: &str) -> Option<Self> {
        match name {
            "chrom" | "CHROM" => Some(Self::Chrom),
            "pos" | "POS" => Some(Self::Pos),
            "id" | "ID" => Some(Self::Id),
            "ref" | "REF" => Some(Self::Ref),
            "alt" | "ALT" => Some(Self::Alt),
            "qual" | "QUAL" => Some(Self::Qual),
            "filter" | "FILTER" => Some(Self::Filter),
            _ => None,
        }
    }

    pub fn arrow_name(self) -> &'static str {
        match self {
            Self::Chrom => "chrom",
            Self::Pos => "pos",
            Self::Id => "id",
            Self::Ref => "ref",
            Self::Alt => "alt",
            Self::Qual => "qual",
            Self::Filter => "filter",
        }
    }
}

pub const DEFAULT_VCF_COLUMNS: [VcfColumn; 7] = [
    VcfColumn::Chrom,
    VcfColumn::Pos,
    VcfColumn::Id,
    VcfColumn::Ref,
    VcfColumn::Alt,
    VcfColumn::Qual,
    VcfColumn::Filter,
];

/// Query options for VCF scanning.
#[derive(Debug, Clone)]
pub struct VcfScanOptions {
    pub columns: Vec<VcfColumn>,
    pub region: Option<FetchRegion>,
}

impl Default for VcfScanOptions {
    fn default() -> Self {
        Self {
            columns: DEFAULT_VCF_COLUMNS.to_vec(),
            region: None,
        }
    }
}

impl VcfScanOptions {
    pub fn wants_column(&self, column: VcfColumn) -> bool {
        self.columns.iter().any(|c| *c == column)
    }
}

/// Columnar VCF table produced by scanners.
#[derive(Debug, Clone, Default)]
pub struct VcfTable {
    pub chrom: Vec<String>,
    pub pos: Vec<i32>,
    pub id: Vec<String>,
    pub reference: Vec<String>,
    pub alt: Vec<String>,
    pub qual: Vec<Option<f32>>,
    pub filter: Vec<String>,
    pub columns: Vec<VcfColumn>,
}

impl VcfTable {
    pub fn new(columns: Vec<VcfColumn>) -> Self {
        Self::with_capacity(columns, 0)
    }

    pub fn with_capacity(columns: Vec<VcfColumn>, rows: usize) -> Self {
        Self {
            columns,
            chrom: Vec::with_capacity(rows),
            pos: Vec::with_capacity(rows),
            id: Vec::with_capacity(rows),
            reference: Vec::with_capacity(rows),
            alt: Vec::with_capacity(rows),
            qual: Vec::with_capacity(rows),
            filter: Vec::with_capacity(rows),
        }
    }

    pub fn len(&self) -> usize {
        self.chrom.len()
    }

    pub fn is_empty(&self) -> bool {
        self.chrom.is_empty()
    }
}