use crate::error::NoodlesError;
use crate::record::{bam_cigar_to_string, bam_quality_to_optional_string, bam_sequence_to_optional_string};
use bamboo_core::{BamColumn, BamScanOptions, BamTable, TagValue};
use bamboo_io::{BamSource, BamStorage};
use noodles::bam as bam;
use noodles::bgzf as bgzf;
use noodles::sam::Header;
use std::fs::File;
use std::io::{BufReader, Cursor, Read};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::Arc;

/// Scan a reader straight into columnar storage without per-record Python objects.
pub fn scan_reader_columnar(
    source: &BamSource,
    header: &Header,
    options: BamScanOptions,
) -> Result<BamTable, NoodlesError> {
    let estimated = estimate_row_capacity(source, options.region.is_some());
    let mut table = BamTable::with_capacity(options.columns.clone(), options.tags.clone(), estimated);

    if options.region.is_some() && bamboo_io::has_index(source) {
        scan_indexed_columnar(source, header, &options, &mut table)?;
    } else {
        scan_linear_columnar(source, header, &options, &mut table)?;
    }

    Ok(table)
}

fn estimate_row_capacity(source: &BamSource, indexed: bool) -> usize {
    if indexed {
        return 16_384;
    }

    match &source.storage {
        BamStorage::Local(path) => path
            .metadata()
            .ok()
            .map(|meta| (meta.len() as usize / 180).clamp(1_024, 8_000_000))
            .unwrap_or(16_384),
        BamStorage::Remote { data, .. } => (data.len() / 180).clamp(1_024, 8_000_000),
    }
}

fn scan_linear_columnar(
    source: &BamSource,
    header: &Header,
    options: &BamScanOptions,
    table: &mut BamTable,
) -> Result<(), NoodlesError> {
    match &source.storage {
        BamStorage::Local(path) => scan_local_linear_columnar(path, header, options, table),
        BamStorage::Remote { data, .. } => scan_remote_linear_columnar(data, header, options, table),
    }
}

fn scan_local_linear_columnar(
    path: &Path,
    header: &Header,
    options: &BamScanOptions,
    table: &mut BamTable,
) -> Result<(), NoodlesError> {
    let workers = std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    let file = File::open(path).map_err(NoodlesError::from)?;
    let decoder = bgzf::MultithreadedReader::with_worker_count(workers, BufReader::new(file));
    let mut reader = bam::io::Reader::from(decoder);
    reader.read_header().map_err(NoodlesError::from)?;
    drain_reader_columnar(&mut reader, header, options, table)
}

fn scan_remote_linear_columnar(
    data: &Arc<[u8]>,
    header: &Header,
    options: &BamScanOptions,
    table: &mut BamTable,
) -> Result<(), NoodlesError> {
    let mut reader = bam::io::Reader::new(Cursor::new(data.clone()));
    reader.read_header().map_err(NoodlesError::from)?;
    drain_reader_columnar(&mut reader, header, options, table)
}

fn scan_indexed_columnar(
    source: &BamSource,
    _header: &Header,
    options: &BamScanOptions,
    table: &mut BamTable,
) -> Result<(), NoodlesError> {
    let region = options
        .region
        .as_ref()
        .ok_or_else(|| NoodlesError::Message("indexed scan requires a region".to_string()))?;
    let parsed_region: noodles::core::Region = region
        .to_samtools_region()
        .parse()
        .map_err(|err: noodles::core::region::ParseError| NoodlesError::Message(err.to_string()))?;

    match &source.storage {
        BamStorage::Local(path) => {
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .map_err(NoodlesError::from)?;
            let header = reader.read_header().map_err(NoodlesError::from)?;
            for result in reader
                .query(&header, &parsed_region)
                .map_err(NoodlesError::from)?
            {
                let record = result.map_err(NoodlesError::from)?;
                if passes_bam_record_filters(&header, &record, options)? {
                    append_bam_record_to_table(&header, &record, table, options);
                }
            }
        }
        BamStorage::Remote { data, index_data, .. } => {
            let index_data = index_data.as_ref().ok_or_else(|| NoodlesError::MissingIndex {
                path: bamboo_io::index_uri_candidates(&source.uri).join(", "),
            })?;
            let mut index_reader = bam::bai::io::Reader::new(index_data.as_ref());
            let index = index_reader.read_index().map_err(NoodlesError::from)?;
            let mut reader = bam::io::indexed_reader::Builder::default()
                .set_index(index)
                .build_from_reader(Cursor::new(data.clone()))
                .map_err(NoodlesError::from)?;
            let header = reader.read_header().map_err(NoodlesError::from)?;
            for result in reader
                .query(&header, &parsed_region)
                .map_err(NoodlesError::from)?
            {
                let record = result.map_err(NoodlesError::from)?;
                if passes_bam_record_filters(&header, &record, options)? {
                    append_bam_record_to_table(&header, &record, table, options);
                }
            }
        }
    }

    Ok(())
}

fn drain_reader_columnar<R: Read>(
    reader: &mut bam::io::Reader<R>,
    header: &Header,
    options: &BamScanOptions,
    table: &mut BamTable,
) -> Result<(), NoodlesError> {
    let mut record = bam::Record::default();
    loop {
        let block_size = reader.read_record(&mut record).map_err(NoodlesError::from)?;
        if block_size == 0 {
            break;
        }
        if passes_bam_record_filters(header, &record, options)? {
            append_bam_record_to_table(header, &record, table, options);
        }
    }
    Ok(())
}

pub fn passes_bam_record_filters(
    header: &Header,
    record: &bam::Record,
    options: &BamScanOptions,
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
        .transpose()
        .map_err(NoodlesError::from)?
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(name, _)| name.to_string());

    if let Some(expected) = &options.reference_name {
        match &reference_name {
            Some(name) if name == expected => {}
            _ => return Ok(false),
        }
    }

    if let Some(region) = &options.region {
        let record_ref = reference_name.as_deref().unwrap_or("");
        if record_ref != region.reference_name {
            return Ok(false);
        }

        let pos = record
            .alignment_start()
            .transpose()
            .map_err(NoodlesError::from)?
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

    Ok(true)
}

pub fn append_bam_record_to_table(
    header: &Header,
    record: &bam::Record,
    table: &mut BamTable,
    options: &BamScanOptions,
) {
    let reference_name = record
        .reference_sequence_id()
        .transpose()
        .ok()
        .flatten()
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(name, _)| name.to_string());

    let mate_reference_name = record
        .mate_reference_sequence_id()
        .transpose()
        .ok()
        .flatten()
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(name, _)| name.to_string());

    let reference_start = record
        .alignment_start()
        .transpose()
        .ok()
        .flatten()
        .map(|position| position.get() as i64 - 1);

    let mate_reference_start = record
        .mate_alignment_start()
        .transpose()
        .ok()
        .flatten()
        .map(|position| position.get() as i64 - 1);

    for column in &options.columns {
        match column {
            BamColumn::QueryName => table.qname.push(record.name().map(|name| name.to_string())),
            BamColumn::Flag => table.flag.push(record.flags().bits()),
            BamColumn::ReferenceName => table.rname.push(reference_name.clone()),
            BamColumn::Position => table.pos.push(reference_start.map(|value| value as i32)),
            BamColumn::MappingQuality => {
                table.mapq.push(record.mapping_quality().map(|quality| quality.get()))
            }
            BamColumn::Cigar => table.cigar.push(bam_cigar_to_string(record.cigar())),
            BamColumn::MateReferenceName => table.rnext.push(mate_reference_name.clone()),
            BamColumn::MatePosition => {
                table.pnext.push(mate_reference_start.map(|value| value as i32))
            }
            BamColumn::TemplateLength => table.tlen.push(Some(record.template_length())),
            BamColumn::Sequence => table.seq.push(bam_sequence_to_optional_string(record.sequence())),
            BamColumn::Quality => {
                table.qual.push(bam_quality_to_optional_string(record.quality_scores()))
            }
        }
    }

    for tag in &mut table.tags {
        tag.values.push(TagValue::Missing);
    }
}