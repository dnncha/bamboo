use crate::options::BamColumn;
use crate::schema::bam_column_name;

/// Typed auxiliary tag values collected during scanning.
#[derive(Debug, Clone, PartialEq)]
pub enum TagValue {
    Int(i64),
    Float(f64),
    String(String),
    Missing,
}

/// One auxiliary tag column.
#[derive(Debug, Clone, PartialEq)]
pub struct TagColumn {
    pub name: String,
    pub values: Vec<TagValue>,
}

/// Columnar BAM table produced by scanners.
#[derive(Debug, Clone, Default)]
pub struct BamTable {
    pub qname: Vec<Option<String>>,
    pub flag: Vec<u16>,
    pub rname: Vec<Option<String>>,
    pub pos: Vec<Option<i32>>,
    pub mapq: Vec<Option<u8>>,
    pub cigar: Vec<String>,
    pub rnext: Vec<Option<String>>,
    pub pnext: Vec<Option<i32>>,
    pub tlen: Vec<Option<i32>>,
    pub seq: Vec<Option<String>>,
    pub qual: Vec<Option<String>>,
    pub tags: Vec<TagColumn>,
    pub columns: Vec<BamColumn>,
}

impl BamTable {
    pub fn new(columns: Vec<BamColumn>, tags: Vec<String>) -> Self {
        Self::with_capacity(columns, tags, 0)
    }

    pub fn with_capacity(columns: Vec<BamColumn>, tags: Vec<String>, rows: usize) -> Self {
        let tag_columns = tags
            .iter()
            .map(|name| TagColumn {
                name: name.clone(),
                values: Vec::with_capacity(rows),
            })
            .collect();

        let mut table = Self {
            columns,
            tags: tag_columns,
            ..Default::default()
        };
        table.reserve(rows);
        table
    }

    pub fn reserve(&mut self, additional: usize) {
        self.qname.reserve(additional);
        self.flag.reserve(additional);
        self.rname.reserve(additional);
        self.pos.reserve(additional);
        self.mapq.reserve(additional);
        self.cigar.reserve(additional);
        self.rnext.reserve(additional);
        self.pnext.reserve(additional);
        self.tlen.reserve(additional);
        self.seq.reserve(additional);
        self.qual.reserve(additional);
        for tag in &mut self.tags {
            tag.values.reserve(additional);
        }
    }

    pub fn len(&self) -> usize {
        let mut rows = 0usize;
        for column in &self.columns {
            rows = rows.max(match column {
                BamColumn::QueryName => self.qname.len(),
                BamColumn::Flag => self.flag.len(),
                BamColumn::ReferenceName => self.rname.len(),
                BamColumn::Position => self.pos.len(),
                BamColumn::MappingQuality => self.mapq.len(),
                BamColumn::Cigar => self.cigar.len(),
                BamColumn::MateReferenceName => self.rnext.len(),
                BamColumn::MatePosition => self.pnext.len(),
                BamColumn::TemplateLength => self.tlen.len(),
                BamColumn::Sequence => self.seq.len(),
                BamColumn::Quality => self.qual.len(),
            });
        }
        for tag in &self.tags {
            rows = rows.max(tag.values.len());
        }
        rows
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn column_names(&self) -> Vec<&'static str> {
        let mut names: Vec<&'static str> = self
            .columns
            .iter()
            .map(|column| bam_column_name(*column))
            .collect();
        for tag in &self.tags {
            names.push(Box::leak(tag.name.clone().into_boxed_str()));
        }
        names
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::options::BamColumn;

    #[test]
    fn len_counts_selected_columns_without_flag() {
        let mut table = BamTable::new(
            vec![
                BamColumn::QueryName,
                BamColumn::ReferenceName,
                BamColumn::Position,
            ],
            Vec::new(),
        );
        table.qname.push(Some("read1".to_string()));
        table.rname.push(Some("chr1".to_string()));
        table.pos.push(Some(99));

        assert_eq!(table.len(), 1);
    }
}