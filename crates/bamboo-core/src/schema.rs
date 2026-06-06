use crate::options::BamColumn;

/// Default fixed columns exported by `scan_bam` / `read_bam`.
pub const DEFAULT_BAM_COLUMNS: [BamColumn; 11] = [
    BamColumn::QueryName,
    BamColumn::Flag,
    BamColumn::ReferenceName,
    BamColumn::Position,
    BamColumn::MappingQuality,
    BamColumn::Cigar,
    BamColumn::MateReferenceName,
    BamColumn::MatePosition,
    BamColumn::TemplateLength,
    BamColumn::Sequence,
    BamColumn::Quality,
];

/// Arrow / Polars column name for a fixed BAM field.
pub fn bam_column_name(column: BamColumn) -> &'static str {
    match column {
        BamColumn::QueryName => "qname",
        BamColumn::Flag => "flag",
        BamColumn::ReferenceName => "rname",
        BamColumn::Position => "pos",
        BamColumn::MappingQuality => "mapq",
        BamColumn::Cigar => "cigar",
        BamColumn::MateReferenceName => "rnext",
        BamColumn::MatePosition => "pnext",
        BamColumn::TemplateLength => "tlen",
        BamColumn::Sequence => "seq",
        BamColumn::Quality => "qual",
    }
}