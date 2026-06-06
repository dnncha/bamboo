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
        let tag_columns = tags
            .into_iter()
            .map(|name| TagColumn {
                name,
                values: Vec::new(),
            })
            .collect();

        Self {
            columns,
            tags: tag_columns,
            ..Default::default()
        }
    }

    pub fn len(&self) -> usize {
        self.flag.len()
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