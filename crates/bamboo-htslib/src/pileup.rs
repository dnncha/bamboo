//! BAM/CRAM pileup via rust-htslib.

use crate::HtslibError;
use rust_htslib::bam::{self, pileup::Pileups, Read};

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

/// Streaming pileup over a fetched region.
pub struct PileupStream {
    #[allow(dead_code)]
    reader: Box<bam::IndexedReader>,
    pileups: Pileups<'static, bam::IndexedReader>,
    target_names: Vec<String>,
    contig: String,
    materialize_reads: bool,
}

impl PileupStream {
    pub fn open(
        path: &str,
        contig: &str,
        start: u32,
        end: u32,
        reference_filename: Option<&str>,
        materialize_reads: bool,
    ) -> Result<Self, HtslibError> {
        let mut reader = Box::new(bam::IndexedReader::from_path(path)?);
        if let Some(reference_path) = reference_filename {
            reader.set_reference(reference_path)?;
        }
        reader.set_threads(4)?;
        reader.fetch((contig, start as i64, end as i64))?;

        let target_names = reader
            .header()
            .target_names()
            .iter()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect::<Vec<_>>();

        let pileups = unsafe {
            let reader_ptr: *mut bam::IndexedReader = reader.as_mut();
            let pileups = (*reader_ptr).pileup();
            std::mem::transmute::<Pileups<'_, bam::IndexedReader>, Pileups<'static, bam::IndexedReader>>(
                pileups,
            )
        };

        Ok(Self {
            reader,
            pileups,
            target_names,
            contig: contig.to_string(),
            materialize_reads,
        })
    }

    pub fn next_column(&mut self) -> Option<Result<PileupColumn, HtslibError>> {
        match self.pileups.next() {
            Some(Ok(pileup)) => Some(Ok(column_from_pileup(
                &pileup,
                &self.target_names,
                &self.contig,
                self.materialize_reads,
            ))),
            Some(Err(err)) => Some(Err(err.into())),
            None => None,
        }
    }
}

fn column_from_pileup(
    pileup: &bam::pileup::Pileup,
    target_names: &[String],
    contig: &str,
    materialize_reads: bool,
) -> PileupColumn {
    let tid = pileup.tid() as i32;
    let reference_name = target_names
        .get(tid as usize)
        .cloned()
        .unwrap_or_else(|| contig.to_string());

    let reads = if materialize_reads {
        pileup
            .alignments()
            .map(|alignment| {
                let record = alignment.record();
                PileupRead {
                    query_name: Some(String::from_utf8_lossy(record.qname()).into_owned()),
                    query_position: alignment.qpos().map(|value| value as u32),
                    is_del: alignment.is_del(),
                    is_head: alignment.is_head(),
                    is_tail: alignment.is_tail(),
                    is_refskip: alignment.is_refskip(),
                }
            })
            .collect()
    } else {
        Vec::new()
    };

    PileupColumn {
        reference_name,
        reference_id: tid,
        position: pileup.pos(),
        depth: pileup.depth(),
        reads,
    }
}

pub fn pileup_region(
    path: &str,
    contig: &str,
    start: u32,
    end: u32,
    reference_filename: Option<&str>,
) -> Result<Vec<PileupColumn>, HtslibError> {
    let mut stream = PileupStream::open(path, contig, start, end, reference_filename, true)?;
    let mut columns = Vec::new();
    while let Some(result) = stream.next_column() {
        columns.push(result?);
    }
    Ok(columns)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bamboo_noodles::fixtures::{
        tiny_bam_path, tiny_cram_path, tiny_fasta_path, write_tiny_bam, write_tiny_bam_index,
        write_tiny_cram, write_tiny_cram_index, write_tiny_fasta,
    };
    use tempfile::tempdir;

    fn assert_read1_pileup(columns: &[PileupColumn]) {
        assert_eq!(columns.len(), 6);
        assert_eq!(columns[0].position, 99);
        assert_eq!(columns[0].depth, 1);
        assert_eq!(columns[0].reads[0].query_name.as_deref(), Some("read1"));
        assert_eq!(columns.last().map(|column| column.position), Some(104));
    }

    fn collect_stream(stream: &mut PileupStream) -> Vec<PileupColumn> {
        let mut columns = Vec::new();
        while let Some(result) = stream.next_column() {
            columns.push(result.unwrap());
        }
        columns
    }

    #[test]
    fn pileup_region_reports_alignment_overlapping_interval() {
        let dir = tempdir().unwrap();
        let bam_path = tiny_bam_path(dir.path());
        write_tiny_bam(&bam_path).unwrap();
        write_tiny_bam_index(&bam_path).unwrap();

        let columns = pileup_region(bam_path.to_str().unwrap(), "chr1", 99, 101, None).unwrap();
        assert_read1_pileup(&columns);
    }

    #[test]
    fn streaming_pileup_matches_materialized_region() {
        let dir = tempdir().unwrap();
        let bam_path = tiny_bam_path(dir.path());
        write_tiny_bam(&bam_path).unwrap();
        write_tiny_bam_index(&bam_path).unwrap();

        let materialized =
            pileup_region(bam_path.to_str().unwrap(), "chr1", 99, 101, None).unwrap();
        let mut stream =
            PileupStream::open(bam_path.to_str().unwrap(), "chr1", 99, 101, None, true).unwrap();
        let streamed = collect_stream(&mut stream);

        assert_eq!(
            streamed
                .iter()
                .map(|column| (column.position, column.depth, column.reads.len()))
                .collect::<Vec<_>>(),
            materialized
                .iter()
                .map(|column| (column.position, column.depth, column.reads.len()))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn streaming_without_reads_reports_depth_only() {
        let dir = tempdir().unwrap();
        let bam_path = tiny_bam_path(dir.path());
        write_tiny_bam(&bam_path).unwrap();
        write_tiny_bam_index(&bam_path).unwrap();

        let mut stream =
            PileupStream::open(bam_path.to_str().unwrap(), "chr1", 99, 101, None, false).unwrap();
        let columns = collect_stream(&mut stream);

        assert_eq!(columns.len(), 6);
        assert!(columns.iter().all(|column| column.reads.is_empty()));
        assert_eq!(columns[0].depth, 1);
    }

    #[test]
    fn cram_pileup_with_external_reference_matches_bam() {
        let dir = tempdir().unwrap();
        let bam_path = tiny_bam_path(dir.path());
        let cram_path = tiny_cram_path(dir.path());
        let fasta_path = tiny_fasta_path(dir.path());
        write_tiny_bam(&bam_path).unwrap();
        write_tiny_bam_index(&bam_path).unwrap();
        write_tiny_cram(&cram_path).unwrap();
        write_tiny_cram_index(&cram_path).unwrap();
        write_tiny_fasta(&fasta_path).unwrap();

        let bam_columns =
            pileup_region(bam_path.to_str().unwrap(), "chr1", 99, 101, None).unwrap();
        let cram_columns = pileup_region(
            cram_path.to_str().unwrap(),
            "chr1",
            99,
            101,
            Some(fasta_path.to_str().unwrap()),
        )
        .unwrap();

        assert_read1_pileup(&cram_columns);
        assert_eq!(
            cram_columns
                .iter()
                .map(|column| (column.position, column.depth))
                .collect::<Vec<_>>(),
            bam_columns
                .iter()
                .map(|column| (column.position, column.depth))
                .collect::<Vec<_>>()
        );
    }
}