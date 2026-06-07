//! BAM pileup via rust-htslib.

use crate::HtslibError;
use rust_htslib::bam::{self, Read};

#[derive(Debug, Clone)]
pub struct PileupColumn {
    pub reference_name: String,
    pub reference_id: i32,
    pub position: u32,
    pub depth: u32,
    pub reads: Vec<PileupRead>,
}

#[derive(Debug, Clone)]
pub struct PileupRead {
    pub query_name: Option<String>,
    pub query_position: Option<u32>,
    pub is_del: bool,
    pub is_head: bool,
    pub is_tail: bool,
    pub is_refskip: bool,
}

pub fn pileup_region(
    path: &str,
    contig: &str,
    start: u32,
    end: u32,
) -> Result<Vec<PileupColumn>, HtslibError> {
    let mut reader = bam::IndexedReader::from_path(path)?;
    reader.fetch((contig, start as i64, end as i64))?;

    let target_names = reader
        .header()
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect::<Vec<_>>();
    let mut columns = Vec::new();

    for pileup in reader.pileup() {
        let pileup = pileup?;
        let tid = pileup.tid() as i32;
        let pos = pileup.pos();

        let reference_name = target_names
            .get(tid as usize)
            .cloned()
            .unwrap_or_else(|| contig.to_string());

        let mut reads = Vec::new();
        for alignment in pileup.alignments() {
            let record = alignment.record();
            reads.push(PileupRead {
                query_name: Some(String::from_utf8_lossy(record.qname()).into_owned()),
                query_position: alignment.qpos().map(|value| value as u32),
                is_del: alignment.is_del(),
                is_head: alignment.is_head(),
                is_tail: alignment.is_tail(),
                is_refskip: alignment.is_refskip(),
            });
        }

        columns.push(PileupColumn {
            reference_name,
            reference_id: tid,
            position: pos,
            depth: pileup.depth(),
            reads,
        });
    }

    Ok(columns)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bamboo_noodles::fixtures::{tiny_bam_path, write_tiny_bam, write_tiny_bam_index};
    use tempfile::tempdir;

    #[test]
    fn pileup_region_reports_alignment_overlapping_interval() {
        let dir = tempdir().unwrap();
        let bam_path = tiny_bam_path(dir.path());
        write_tiny_bam(&bam_path).unwrap();
        write_tiny_bam_index(&bam_path).unwrap();

        let columns = pileup_region(bam_path.to_str().unwrap(), "chr1", 99, 101).unwrap();
        assert_eq!(columns.len(), 6);
        assert_eq!(columns[0].position, 99);
        assert_eq!(columns[0].depth, 1);
        assert_eq!(columns[0].reads[0].query_name.as_deref(), Some("read1"));
        assert_eq!(columns.last().map(|column| column.position), Some(104));
    }
}