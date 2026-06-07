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

/// Build a tiny plain VCF payload used in tests.
pub fn tiny_vcf_bytes() -> io::Result<Vec<u8>> {
    Ok(
        b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=1000>\n\
##contig=<ID=chr2,length=1000>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t60\tPASS\t.\n\
chr2\t250\trs1\tC\tT\t30\tPASS\tAF=0.5\n"
            .to_vec(),
    )
}

/// Write the tiny VCF fixture to `path`.
pub fn write_tiny_vcf(path: &Path) -> io::Result<()> {
    std::fs::write(path, tiny_vcf_bytes()?)
}

/// Convenience path helper for temp directories.
pub fn tiny_vcf_path(dir: &Path) -> PathBuf {
    dir.join("tiny.vcf")
}

/// Path helper for a bgzipped tiny VCF fixture.
pub fn tiny_vcf_gz_path(dir: &Path) -> PathBuf {
    dir.join("tiny.vcf.gz")
}

/// Path helper for a tiny BCF fixture.
pub fn tiny_bcf_path(dir: &Path) -> PathBuf {
    dir.join("tiny.bcf")
}

/// Path helper for a tiny CRAM fixture.
pub fn tiny_cram_path(dir: &Path) -> PathBuf {
    dir.join("tiny.cram")
}

/// Path helper for the CRAI sidecar of `tiny.cram`.
pub fn tiny_cram_index_path(dir: &Path) -> PathBuf {
    let mut path = tiny_cram_path(dir).into_os_string();
    path.push(".crai");
    PathBuf::from(path)
}

/// Path helper for a tiny FASTA reference fixture.
pub fn tiny_fasta_path(dir: &Path) -> PathBuf {
    dir.join("tiny.fasta")
}

/// Write a tiny CRAM fixture with two mapped reads.
pub fn write_tiny_cram(path: &Path) -> io::Result<()> {
    use noodles::cram as cram;
    use noodles::fasta as fasta;
    use noodles::fasta::record::{Definition, Sequence as FastaSequence};
    use noodles::sam::alignment::io::Write as _;

    let reference_sequences = vec![
        fasta::Record::new(
            Definition::new("chr1", None),
            FastaSequence::from(vec![b'N'; 1000]),
        ),
        fasta::Record::new(
            Definition::new("chr2", None),
            FastaSequence::from(vec![b'N'; 1000]),
        ),
    ];

    let mut header_builder = Header::builder().set_header(
        Map::<HeaderMap>::builder()
            .insert(SORT_ORDER, COORDINATE)
            .build()
            .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.to_string()))?,
    );

    for record in &reference_sequences {
        let length = NonZeroUsize::new(record.sequence().len()).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "invalid reference length")
        })?;
        header_builder = header_builder.add_reference_sequence(
            record.name(),
            Map::<ReferenceSequence>::new(length),
        );
    }

    let header = header_builder.build();
    let repository = fasta::Repository::new(reference_sequences);
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(path)?;
    writer.write_header(&header)?;

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(noodles::core::Position::try_from(100).unwrap())
        .set_mapping_quality(MappingQuality::try_from(60).unwrap())
        .set_cigar([Op::new(Kind::Match, 6)].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(b"ACGTAC".to_vec()))
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
        .set_sequence(Sequence::from(b"TTTT".to_vec()))
        .set_quality_scores(QualityScores::from(Vec::from(b"!!!!".as_slice())))
        .build();
    writer.write_alignment_record(&header, &record2)?;

    writer.try_finish(&header)?;
    Ok(())
}

/// Write a CRAI index for a tiny CRAM fixture.
pub fn write_tiny_cram_index(cram_path: &Path) -> io::Result<()> {
    use noodles::cram::{self as cram, crai};

    let index = cram::fs::index(cram_path)?;
    let mut index_path = cram_path.as_os_str().to_os_string();
    index_path.push(".crai");

    let file = std::fs::File::create(PathBuf::from(index_path))?;
    let mut writer = crai::io::Writer::new(file);
    writer.write_index(&index)?;
    writer.finish()?;
    Ok(())
}

/// Write a tiny FASTA reference fixture matching the CRAM/BAM headers.
pub fn write_tiny_fasta(path: &Path) -> io::Result<()> {
    use noodles::fasta as fasta;
    use noodles::fasta::record::{Definition, Sequence};

    let records = [
        fasta::Record::new(
            Definition::new("chr1", None),
            Sequence::from(vec![b'N'; 1000]),
        ),
        fasta::Record::new(
            Definition::new("chr2", None),
            Sequence::from(vec![b'N'; 1000]),
        ),
    ];

    let mut writer = fasta::io::writer::Builder::default().build_from_path(path)?;
    for record in records {
        writer.write_record(&record)?;
    }
    Ok(())
}

/// Write a CSI index for a tiny BCF fixture.
pub fn write_tiny_bcf_index(bcf_path: &Path) -> io::Result<()> {
    crate::bcf::write_bcf_index(bcf_path)
}

/// Write a tiny BCF fixture equivalent to `tiny.vcf`.
pub fn write_tiny_bcf(path: &Path) -> io::Result<()> {
    use noodles::bcf as bcf;
    use noodles::vcf as vcf;
    use noodles::vcf::variant::io::Write as _;

    let vcf_path = path
        .parent()
        .map(|dir| dir.join("tiny.vcf"))
        .filter(|candidate| candidate.exists())
        .unwrap_or_else(|| path.with_extension("vcf"));

    if !vcf_path.exists() {
        write_tiny_vcf(&vcf_path)?;
    }

    let mut vcf_reader = vcf::io::reader::Builder::default().build_from_path(&vcf_path)?;
    let header = vcf_reader.read_header()?;

    let mut bcf_writer = bcf::io::writer::Builder::default().build_from_path(path)?;
    bcf_writer.write_header(&header)?;

    for result in vcf_reader.record_bufs(&header) {
        let record = result?;
        bcf_writer.write_variant_record(&header, &record)?;
    }

    Ok(())
}

/// Write a bgzipped tiny VCF fixture to `path`.
pub fn write_tiny_vcf_gz(path: &Path) -> io::Result<()> {
    use std::io::Write as _;
    let file = std::fs::File::create(path)?;
    let mut writer = noodles::bgzf::Writer::new(file);
    writer.write_all(&tiny_vcf_bytes()?)?;
    writer.try_finish()?;
    Ok(())
}

/// Write a tabix index for a bgzipped VCF fixture.
pub fn write_tiny_vcf_index(vcf_gz_path: &Path) -> io::Result<()> {
    let index = noodles::vcf::fs::index(vcf_gz_path)?;
    let mut index_path = std::ffi::OsString::from(vcf_gz_path);
    index_path.push(".tbi");
    noodles::tabix::fs::write(PathBuf::from(index_path), &index)
}

const BENCH_REFERENCES: [(&str, usize); 5] = [
    ("chr1", 250_000_000),
    ("chr2", 242_193_529),
    ("chr3", 198_295_559),
    ("chr4", 190_214_555),
    ("chr5", 181_538_259),
];

const BENCH_READ_LENGTH: usize = 100;

fn coordinate_sorted_header(
    references: &[(&str, usize)],
) -> io::Result<Header> {
    let mut builder = Header::builder().set_header(
        Map::<HeaderMap>::builder()
            .insert(SORT_ORDER, COORDINATE)
            .build()
            .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.to_string()))?,
    );

    for (name, length) in references {
        let length = NonZeroUsize::new(*length).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("reference '{name}' has invalid length {length}"),
            )
        })?;
        builder = builder.add_reference_sequence(
            *name,
            Map::<ReferenceSequence>::new(length),
        );
    }

    Ok(builder.build())
}

/// Write a coordinate-sorted synthetic BAM for benchmarks.
pub fn write_bench_bam(path: &Path, record_count: usize) -> io::Result<()> {
    let mut writer = bam::io::writer::Builder::default().build_from_path(path)?;
    let header = coordinate_sorted_header(&BENCH_REFERENCES)?;
    writer.write_header(&header)?;

    let sequence = Sequence::from(vec![b'A'; BENCH_READ_LENGTH]);
    let quality_scores = QualityScores::from(vec![b'!'; BENCH_READ_LENGTH]);
    let cigar: Cigar = [Op::new(Kind::Match, BENCH_READ_LENGTH)]
        .into_iter()
        .collect();

    let refs = BENCH_REFERENCES.len();
    let base_per_ref = record_count / refs;
    let remainder = record_count % refs;
    let mut global_index = 0usize;

    for reference_sequence_id in 0..refs {
        let records_for_ref = base_per_ref + usize::from(reference_sequence_id < remainder);
        for offset in 0..records_for_ref {
            let alignment_start = 1_000 + offset * 250;
            let mapping_quality = ((global_index % 60) + 1) as u8;
            let name = format!("bench_read_{global_index:08}");

            let record = RecordBuf::builder()
                .set_name(name.as_str())
                .set_flags(Flags::empty())
                .set_reference_sequence_id(reference_sequence_id)
                .set_alignment_start(
                    noodles::core::Position::try_from(alignment_start)
                        .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.to_string()))?,
                )
                .set_mapping_quality(
                    MappingQuality::try_from(mapping_quality)
                        .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.to_string()))?,
                )
                .set_cigar(cigar.clone())
                .set_sequence(sequence.clone())
                .set_quality_scores(quality_scores.clone())
                .build();

            writer.write_alignment_record(&header, &record)?;
            global_index += 1;
        }
    }

    writer.try_finish()?;
    Ok(())
}

/// Write a BAI index for a benchmark BAM.
pub fn write_bench_bam_index(bam_path: &Path) -> io::Result<()> {
    let index = bam::fs::index(bam_path)?;
    let index_path = bam_path.with_extension("bam.bai");
    bam::bai::fs::write(index_path, &index)
}

const BENCH_CRAM_REFERENCES: [(&str, usize); 5] = [
    ("chr1", 10_000_000),
    ("chr2", 10_000_000),
    ("chr3", 10_000_000),
    ("chr4", 10_000_000),
    ("chr5", 10_000_000),
];

/// Write a coordinate-sorted synthetic CRAM for benchmarks.
pub fn write_bench_cram(path: &Path, record_count: usize) -> io::Result<()> {
    use noodles::cram as cram;
    use noodles::fasta as fasta;
    use noodles::fasta::record::{Definition, Sequence as FastaSequence};
    use noodles::sam::alignment::io::Write as _;

    let reference_sequences = BENCH_CRAM_REFERENCES
        .iter()
        .map(|(name, length)| {
            fasta::Record::new(
                Definition::new(*name, None),
                FastaSequence::from(vec![b'N'; *length]),
            )
        })
        .collect::<Vec<_>>();

    let header = coordinate_sorted_header(&BENCH_CRAM_REFERENCES)?;
    let repository = fasta::Repository::new(reference_sequences);
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(path)?;
    writer.write_header(&header)?;

    let sequence = Sequence::from(vec![b'A'; BENCH_READ_LENGTH]);
    let quality_scores = QualityScores::from(vec![b'!'; BENCH_READ_LENGTH]);
    let cigar: Cigar = [Op::new(Kind::Match, BENCH_READ_LENGTH)]
        .into_iter()
        .collect();

    let refs = BENCH_CRAM_REFERENCES.len();
    let base_per_ref = record_count / refs;
    let remainder = record_count % refs;
    let mut global_index = 0usize;

    for reference_sequence_id in 0..refs {
        let records_for_ref = base_per_ref + usize::from(reference_sequence_id < remainder);
        for offset in 0..records_for_ref {
            let alignment_start = 1_000 + offset * 250;
            let mapping_quality = ((global_index % 60) + 1) as u8;
            let name = format!("bench_read_{global_index:08}");

            let record = RecordBuf::builder()
                .set_name(name.as_str())
                .set_flags(Flags::empty())
                .set_reference_sequence_id(reference_sequence_id)
                .set_alignment_start(
                    noodles::core::Position::try_from(alignment_start)
                        .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.to_string()))?,
                )
                .set_mapping_quality(
                    MappingQuality::try_from(mapping_quality)
                        .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.to_string()))?,
                )
                .set_cigar(cigar.clone())
                .set_sequence(sequence.clone())
                .set_quality_scores(quality_scores.clone())
                .build();

            writer.write_alignment_record(&header, &record)?;
            global_index += 1;
        }
    }

    writer.try_finish(&header)?;
    Ok(())
}

/// Write a CRAI index for a benchmark CRAM.
pub fn write_bench_cram_index(cram_path: &Path) -> io::Result<()> {
    write_tiny_cram_index(cram_path)
}

/// Write a FAI index for a FASTA fixture.
pub fn write_fasta_index(fasta_path: &Path) -> io::Result<()> {
    use noodles::fasta as fasta;

    let index = fasta::fs::index(fasta_path)?;
    let mut index_path = std::ffi::OsString::from(fasta_path);
    index_path.push(".fai");
    fasta::fai::fs::write(PathBuf::from(index_path), &index)
}

/// Write a FAI sidecar for `tiny.fasta`.
pub fn write_tiny_fasta_index(fasta_path: &Path) -> io::Result<()> {
    write_fasta_index(fasta_path)
}

/// Write a FAI sidecar for a benchmark FASTA.
pub fn write_bench_fasta_index(fasta_path: &Path) -> io::Result<()> {
    write_fasta_index(fasta_path)
}

/// Write a FASTA reference for benchmark CRAMs.
pub fn write_bench_fasta(path: &Path) -> io::Result<()> {
    use noodles::fasta as fasta;
    use noodles::fasta::record::{Definition, Sequence};

    let mut writer = fasta::io::writer::Builder::default().build_from_path(path)?;
    for (name, length) in BENCH_CRAM_REFERENCES {
        writer.write_record(&fasta::Record::new(
            Definition::new(name, None),
            Sequence::from(vec![b'N'; length]),
        ))?;
    }
    Ok(())
}