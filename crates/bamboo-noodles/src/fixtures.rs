//! Test fixtures for Bamboo BAM readers.

use noodles::bam as bam;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER};
use noodles::sam::header::record::value::map::Header as HeaderMap;
use std::io;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

/// Build the shared coordinate-sorted tiny BAM payload used in tests.
pub fn tiny_bam_bytes() -> io::Result<Vec<u8>> {
    let mut writer = bam::io::Writer::new(Vec::new());

    let header = Header::builder()
        .set_header(
            Map::<HeaderMap>::builder()
                .insert(SORT_ORDER, COORDINATE)
                .build()
                .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.to_string()))?,
        )
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
        )
        .add_reference_sequence(
            "chr2",
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
        )
        .build();

    writer.write_header(&header)?;

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(noodles::core::Position::try_from(100).unwrap())
        .set_mapping_quality(MappingQuality::try_from(60).unwrap())
        .set_cigar([Op::new(Kind::Match, 6)].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(b"ACGTAC"))
        .set_quality_scores(QualityScores::from(Vec::from(b"!!!!!!".as_slice())))
        .build();
    writer.write_alignment_record(&header, &record)?;

    let record2 = RecordBuf::builder()
        .set_name("read2")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(1)
        .set_alignment_start(noodles::core::Position::try_from(250).unwrap())
        .set_mapping_quality(MappingQuality::try_from(5).unwrap())
        .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(b"TTTT"))
        .set_quality_scores(QualityScores::from(Vec::from(b"!!!!".as_slice())))
        .build();
    writer.write_alignment_record(&header, &record2)?;

    writer.try_finish()?;
    Ok(writer.into_inner().into_inner())
}

/// Write the tiny coordinate-sorted BAM fixture to `path`.
pub fn write_tiny_bam(path: &Path) -> io::Result<()> {
    std::fs::write(path, tiny_bam_bytes()?)
}

/// Write a BAI index for the tiny BAM fixture.
pub fn write_tiny_bam_index(bam_path: &Path) -> io::Result<()> {
    let index = bam::fs::index(bam_path)?;
    let index_path = bam_path.with_extension("bam.bai");
    bam::bai::fs::write(index_path, &index)
}

/// Convenience path helper for temp directories.
pub fn tiny_bam_path(dir: &Path) -> PathBuf {
    dir.join("tiny.bam")
}