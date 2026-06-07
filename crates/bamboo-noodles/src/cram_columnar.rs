use crate::cram::{cram_index_path, read_header_from_path, reference_repository};
use crate::error::NoodlesError;
use bamboo_core::{BamColumn, BamScanOptions, BamTable, TagValue};
use noodles::cram as cram;
use noodles::sam::Header;
use noodles::sam::alignment::Record as SamAlignmentRecord;
use noodles::sam::alignment::record::cigar::op::Kind;

use std::path::Path;

/// Scan a CRAM file directly into columnar storage without per-record Python objects.
pub fn scan_cram_columnar(
    path: &str,
    options: BamScanOptions,
    reference_fasta: Option<&str>,
) -> Result<BamTable, NoodlesError> {
    let header = read_header_from_path(path)?;
    let repository = reference_repository(&header, reference_fasta)?;
    let estimated = estimate_row_capacity(path, options.region.is_some());
    let mut table = BamTable::with_capacity(options.columns.clone(), options.tags.clone(), estimated);

    if options.region.is_some() && cram_index_path(path).exists() {
        scan_indexed_cram_columnar(path, &repository, &options, &mut table)?;
    } else {
        scan_linear_cram_columnar(path, &repository, &options, &mut table)?;
    }

    Ok(table)
}

fn estimate_row_capacity(path: &str, indexed: bool) -> usize {
    if indexed {
        return 16_384;
    }

    Path::new(path)
        .metadata()
        .ok()
        .map(|meta| (meta.len() as usize / 180).clamp(1_024, 8_000_000))
        .unwrap_or(16_384)
}

fn scan_linear_cram_columnar(
    path: &str,
    repository: &noodles::fasta::Repository,
    options: &BamScanOptions,
    table: &mut BamTable,
) -> Result<(), NoodlesError> {
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
    let header = reader.read_header().map_err(NoodlesError::from)?;

    for result in reader.records(&header) {
        let record = result.map_err(NoodlesError::from)?;
        if passes_cram_record_filters(&header, &record, options, false)? {
            append_cram_record_to_table(&header, &record, table, options);
        }
    }

    Ok(())
}

fn scan_indexed_cram_columnar(
    path: &str,
    repository: &noodles::fasta::Repository,
    options: &BamScanOptions,
    table: &mut BamTable,
) -> Result<(), NoodlesError> {
    let region = options
        .region
        .as_ref()
        .ok_or_else(|| NoodlesError::Message("indexed CRAM scan requires a region".to_string()))?;
    let parsed_region: noodles::core::Region = region
        .to_samtools_region()
        .parse()
        .map_err(|err: noodles::core::region::ParseError| NoodlesError::Message(err.to_string()))?;

    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
    let header = reader.read_header().map_err(NoodlesError::from)?;

    for result in reader
        .query(&header, &parsed_region)
        .map_err(NoodlesError::from)?
    {
        let record = result.map_err(NoodlesError::from)?;
        if passes_cram_record_filters(&header, &record, options, true)? {
            append_cram_record_to_table(&header, &record, table, options);
        }
    }

    Ok(())
}

fn passes_cram_record_filters(
    header: &Header,
    record: &cram::Record,
    options: &BamScanOptions,
    skip_region_filter: bool,
) -> Result<bool, NoodlesError> {
    if let Some(min_mapq) = options.min_mapq {
        let mapq = record.mapping_quality().map(|quality| quality.get());
        match mapq {
            Some(value) if value >= min_mapq => {}
            _ => return Ok(false),
        }
    }

    let reference_name = record
        .reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(name, _)| name.to_string());

    if let Some(expected) = &options.reference_name {
        match &reference_name {
            Some(name) if name == expected => {}
            _ => return Ok(false),
        }
    }

    if !skip_region_filter {
        if let Some(region) = &options.region {
            let record_ref = reference_name.as_deref().unwrap_or("");
            if record_ref != region.reference_name {
                return Ok(false);
            }

            let pos = record
                .alignment_start()
                .map(|position| position.get() as i64 - 1)
                .unwrap_or(-1);

            if let Some(start) = region.start {
                if pos < start as i64 {
                    return Ok(false);
                }
            }
            if let Some(end) = region.end {
                if pos >= end as i64 {
                    return Ok(false);
                }
            }
        }
    }

    Ok(true)
}

fn append_cram_record_to_table(
    header: &Header,
    record: &cram::Record,
    table: &mut BamTable,
    options: &BamScanOptions,
) {
    let reference_name = record
        .reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(name, _)| name.to_string());

    let mate_reference_name = record
        .next_fragment_reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(name, _)| name.to_string());

    let reference_start = record
        .alignment_start()
        .map(|position| position.get() as i64 - 1);

    let mate_reference_start = record
        .next_mate_alignment_start()
        .map(|position| position.get() as i64 - 1);

    for column in &options.columns {
        match column {
            BamColumn::QueryName => {
                table.qname.push(record.name().map(|name| name.to_string()))
            }
            BamColumn::Flag => table.flag.push(record.flags().bits()),
            BamColumn::ReferenceName => table.rname.push(reference_name.clone()),
            BamColumn::Position => table.pos.push(reference_start.map(|value| value as i32)),
            BamColumn::MappingQuality => {
                table.mapq.push(record.mapping_quality().map(|quality| quality.get()))
            }
            BamColumn::Cigar => table.cigar.push(sam_cigar_to_string(record)),
            BamColumn::MateReferenceName => table.rnext.push(mate_reference_name.clone()),
            BamColumn::MatePosition => {
                table.pnext.push(mate_reference_start.map(|value| value as i32))
            }
            BamColumn::TemplateLength => table.tlen.push(Some(record.template_length())),
            BamColumn::Sequence => table.seq.push(sam_sequence_to_optional_string(record)),
            BamColumn::Quality => {
                table.qual.push(sam_quality_to_optional_string(record))
            }
        }
    }

    for tag in &mut table.tags {
        tag.values.push(TagValue::Missing);
    }
}

fn sam_cigar_to_string(record: &cram::Record) -> String {
    SamAlignmentRecord::cigar(record)
        .iter()
        .filter_map(Result::ok)
        .map(|op| format!("{}{}", op.len(), kind_to_char(op.kind())))
        .collect()
}

fn sam_sequence_to_optional_string(record: &cram::Record) -> Option<String> {
    let bytes = record.sequence().as_ref();
    if bytes.is_empty() {
        None
    } else {
        Some(String::from_utf8_lossy(bytes).into_owned())
    }
}

fn sam_quality_to_optional_string(record: &cram::Record) -> Option<String> {
    let bytes = record.quality_scores().as_ref();
    if bytes.is_empty() {
        None
    } else {
        Some(String::from_utf8_lossy(bytes).into_owned())
    }
}

fn kind_to_char(kind: Kind) -> char {
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{
        tiny_bam_path, tiny_cram_path, tiny_fasta_path, write_tiny_bam, write_tiny_cram,
        write_tiny_cram_index, write_tiny_fasta,
    };
    use bamboo_core::{BamColumn, FetchRegion};
    use tempfile::tempdir;

    #[test]
    fn indexed_cram_columnar_matches_linear_region_filter() {
        let dir = tempdir().unwrap();
        let linear_path = dir.path().join("linear.cram");
        let indexed_path = tiny_cram_path(dir.path());
        write_tiny_cram(&linear_path).unwrap();
        write_tiny_cram(&indexed_path).unwrap();
        write_tiny_cram_index(&indexed_path).unwrap();

        let options = BamScanOptions {
            region: Some(FetchRegion::from_samtools_region("chr1:100-101").expect("valid region")),
            columns: vec![
                BamColumn::QueryName,
                BamColumn::ReferenceName,
                BamColumn::Position,
                BamColumn::MappingQuality,
                BamColumn::Cigar,
            ],
            ..Default::default()
        };

        assert!(cram_index_path(indexed_path.to_str().unwrap()).exists());

        let indexed = scan_cram_columnar(indexed_path.to_str().unwrap(), options.clone(), None).unwrap();
        let linear = scan_cram_columnar(linear_path.to_str().unwrap(), options, None).unwrap();

        assert_eq!(indexed.len(), 1);
        assert_eq!(indexed.qname, linear.qname);
        assert_eq!(indexed.rname, linear.rname);
        assert_eq!(indexed.pos, linear.pos);
        assert_eq!(indexed.mapq, linear.mapq);
        assert_eq!(indexed.cigar, linear.cigar);
    }

    #[test]
    fn cram_columnar_with_external_fasta_matches_bam() {
        let dir = tempdir().unwrap();
        let cram_path = tiny_cram_path(dir.path());
        let bam_path = tiny_bam_path(dir.path());
        let fasta_path = tiny_fasta_path(dir.path());
        write_tiny_cram(&cram_path).unwrap();
        write_tiny_bam(&bam_path).unwrap();
        write_tiny_fasta(&fasta_path).unwrap();

        let options = BamScanOptions {
            columns: vec![
                BamColumn::QueryName,
                BamColumn::Position,
                BamColumn::MappingQuality,
                BamColumn::Cigar,
            ],
            ..Default::default()
        };

        let cram_table = scan_cram_columnar(
            cram_path.to_str().unwrap(),
            options.clone(),
            Some(fasta_path.to_str().unwrap()),
        )
        .unwrap();
        let bam_reader = crate::BamReader::open(bam_path.to_str().unwrap()).unwrap();
        let bam_table = crate::scan::scan_reader(&bam_reader, options).unwrap();

        assert_eq!(cram_table.len(), bam_table.len());
        assert_eq!(cram_table.qname, bam_table.qname);
        assert_eq!(cram_table.pos, bam_table.pos);
        assert_eq!(cram_table.mapq, bam_table.mapq);
        assert_eq!(cram_table.cigar, bam_table.cigar);
    }
}