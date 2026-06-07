use crate::error::NoodlesError;
use bamboo_core::{FetchRegion, VcfColumn, VcfScanOptions, VcfTable};
use noodles::vcf as vcf;
use noodles::vcf::variant::record::{AlternateBases, Filters, Ids};
use std::ffi::OsString;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// High-level VCF reader backed by noodles.
pub struct VcfReader {
    path: String,
    header: vcf::Header,
}

impl VcfReader {
    pub fn open(path: &str) -> Result<Self, NoodlesError> {
        let header = read_header_from_path(path)?;
        Ok(Self {
            path: path.to_string(),
            header,
        })
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    pub fn header(&self) -> &vcf::Header {
        &self.header
    }

    pub fn has_index(&self) -> bool {
        has_tabix_index(&self.path)
    }

    pub fn reference_names(&self) -> Vec<String> {
        self.header
            .contigs()
            .keys()
            .map(|name| name.to_string())
            .collect()
    }

    pub fn count_records(&self) -> Result<usize, NoodlesError> {
        scan_vcf(&self.path, VcfScanOptions::default()).map(|table| table.len())
    }

    pub fn scan(&self, options: VcfScanOptions) -> Result<VcfTable, NoodlesError> {
        scan_vcf(&self.path, options)
    }
}

pub fn scan_vcf(path: &str, options: VcfScanOptions) -> Result<VcfTable, NoodlesError> {
    if options.region.is_some() && has_tabix_index(path) {
        scan_vcf_indexed(path, &options)
    } else {
        scan_vcf_linear(path, &options)
    }
}

fn scan_vcf_linear(path: &str, options: &VcfScanOptions) -> Result<VcfTable, NoodlesError> {
    let mut reader = open_reader(path)?;
    let header = reader.read_header().map_err(NoodlesError::from)?;
    let mut table = VcfTable::new(options.columns.clone());

    for result in reader.record_bufs(&header) {
        let record = result.map_err(NoodlesError::from)?;
        if !passes_region_filter_buf(&record, &options.region) {
            continue;
        }
        append_record_buf(&mut table, &record, options);
    }

    Ok(table)
}

fn scan_vcf_indexed(path: &str, options: &VcfScanOptions) -> Result<VcfTable, NoodlesError> {
    let region = options
        .region
        .as_ref()
        .ok_or_else(|| NoodlesError::Message("indexed VCF scan requires a region".to_string()))?;
    let parsed_region: noodles::core::Region = region
        .to_samtools_region()
        .parse()
        .map_err(|err: noodles::core::region::ParseError| NoodlesError::Message(err.to_string()))?;

    let mut reader = vcf::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
    let header = reader.read_header().map_err(NoodlesError::from)?;
    let mut table = VcfTable::new(options.columns.clone());

    for result in reader
        .query(&header, &parsed_region)
        .map_err(NoodlesError::from)?
    {
        let record = result.map_err(NoodlesError::from)?;
        append_vcf_record(&mut table, &header, &record, options)?;
    }

    Ok(table)
}

pub fn has_tabix_index(path: &str) -> bool {
    let path = Path::new(path);
    index_path(path, "tbi").exists() || index_path(path, "csi").exists()
}

fn index_path(path: &Path, ext: &str) -> PathBuf {
    let mut candidate = OsString::from(path);
    candidate.push(".");
    candidate.push(ext);
    PathBuf::from(candidate)
}

fn read_header_from_path(path: &str) -> Result<vcf::Header, NoodlesError> {
    let mut reader = open_reader(path)?;
    reader.read_header().map_err(NoodlesError::from)
}

fn open_reader(path: &str) -> Result<vcf::io::Reader<Box<dyn BufRead>>, NoodlesError> {
    vcf::io::reader::Builder::default()
        .build_from_path(Path::new(path))
        .map_err(NoodlesError::from)
}

pub(crate) fn passes_region_filter_buf(
    record: &vcf::variant::RecordBuf,
    region: &Option<FetchRegion>,
) -> bool {
    let Some(region) = region else {
        return true;
    };

    if record.reference_sequence_name() != region.reference_name {
        return false;
    }

    let pos = record
        .variant_start()
        .map(|position| position.get() as i64)
        .unwrap_or(-1);

    region_contains_position(region, pos)
}

fn region_contains_position(region: &FetchRegion, pos: i64) -> bool {
    if let Some(start) = region.start {
        if pos < start as i64 {
            return false;
        }
    }
    if let Some(end) = region.end {
        if pos >= end as i64 {
            return false;
        }
    }
    true
}

pub(crate) fn append_record_buf(
    table: &mut VcfTable,
    record: &vcf::variant::RecordBuf,
    options: &VcfScanOptions,
) {
    let fields = extract_record_buf_fields(record);
    append_fields(table, &fields, options);
}

pub(crate) fn append_vcf_record(
    table: &mut VcfTable,
    header: &vcf::Header,
    record: &vcf::Record,
    options: &VcfScanOptions,
) -> Result<(), NoodlesError> {
    let fields = extract_vcf_record_fields(header, record)?;
    append_fields(table, &fields, options);
    Ok(())
}

struct VcfRecordFields {
    chrom: String,
    pos: i32,
    id: String,
    reference: String,
    alt: String,
    qual: Option<f32>,
    filter: String,
}

fn extract_record_buf_fields(record: &vcf::variant::RecordBuf) -> VcfRecordFields {
    VcfRecordFields {
        chrom: record.reference_sequence_name().to_string(),
        pos: record
            .variant_start()
            .map(|position| position.get() as i32)
            .unwrap_or(0),
        id: record
            .ids()
            .iter()
            .next()
            .unwrap_or(".")
            .to_string(),
        reference: record.reference_bases().to_string(),
        alt: record
            .alternate_bases()
            .as_ref()
            .iter()
            .map(|allele| allele.as_str())
            .collect::<Vec<_>>()
            .join(","),
        qual: record.quality_score(),
        filter: record
            .filters()
            .as_ref()
            .iter()
            .next()
            .cloned()
            .unwrap_or_else(|| ".".to_string()),
    }
}

fn extract_vcf_record_fields(
    header: &vcf::Header,
    record: &vcf::Record,
) -> Result<VcfRecordFields, NoodlesError> {
    let pos = match record.variant_start() {
        Some(Ok(position)) => position.get() as i32,
        Some(Err(err)) => return Err(NoodlesError::from(err)),
        None => 0,
    };
    let qual = match record.quality_score() {
        Some(Ok(score)) => Some(score),
        Some(Err(err)) => return Err(NoodlesError::from(err)),
        None => None,
    };
    let id = record.ids().iter().next().unwrap_or(".").to_string();
    let alt = record
        .alternate_bases()
        .iter()
        .map(|result| result.map(|allele| allele.to_string()))
        .collect::<Result<Vec<_>, _>>()
        .map_err(NoodlesError::from)?
        .join(",");
    let filter = record
        .filters()
        .iter(header)
        .next()
        .transpose()
        .map_err(NoodlesError::from)?
        .unwrap_or(".")
        .to_string();

    Ok(VcfRecordFields {
        chrom: record.reference_sequence_name().to_string(),
        pos,
        id,
        reference: record.reference_bases().to_string(),
        alt,
        qual,
        filter,
    })
}

fn append_fields(table: &mut VcfTable, fields: &VcfRecordFields, options: &VcfScanOptions) {
    if options.wants_column(VcfColumn::Chrom) {
        table.chrom.push(fields.chrom.clone());
    }
    if options.wants_column(VcfColumn::Pos) {
        table.pos.push(fields.pos);
    }
    if options.wants_column(VcfColumn::Id) {
        table.id.push(fields.id.clone());
    }
    if options.wants_column(VcfColumn::Ref) {
        table.reference.push(fields.reference.clone());
    }
    if options.wants_column(VcfColumn::Alt) {
        table.alt.push(fields.alt.clone());
    }
    if options.wants_column(VcfColumn::Qual) {
        table.qual.push(fields.qual);
    }
    if options.wants_column(VcfColumn::Filter) {
        table.filter.push(fields.filter.clone());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{
        tiny_vcf_gz_path, tiny_vcf_path, write_tiny_vcf, write_tiny_vcf_gz, write_tiny_vcf_index,
    };
    use tempfile::tempdir;

    #[test]
    fn reads_tiny_vcf() {
        let dir = tempdir().unwrap();
        let path = tiny_vcf_path(dir.path());
        write_tiny_vcf(&path).unwrap();

        let reader = VcfReader::open(path.to_str().unwrap()).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
        assert_eq!(reader.reference_names(), vec!["chr1", "chr2"]);
    }

    #[test]
    fn filters_by_region() {
        let dir = tempdir().unwrap();
        let path = tiny_vcf_path(dir.path());
        write_tiny_vcf(&path).unwrap();

        let options = VcfScanOptions {
            region: Some(FetchRegion {
                reference_name: "chr1".to_string(),
                start: Some(99),
                end: Some(200),
            }),
            ..Default::default()
        };
        let table = scan_vcf(path.to_str().unwrap(), options).unwrap();
        assert_eq!(table.len(), 1);
        assert_eq!(table.chrom[0], "chr1");
        assert_eq!(table.pos[0], 100);
    }

    #[test]
    fn indexed_fetch_matches_linear_scan() {
        let dir = tempdir().unwrap();
        let path = tiny_vcf_gz_path(dir.path());
        write_tiny_vcf_gz(&path).unwrap();
        write_tiny_vcf_index(&path).unwrap();

        let region = VcfScanOptions {
            region: Some(FetchRegion {
                reference_name: "chr1".to_string(),
                start: Some(99),
                end: Some(200),
            }),
            ..Default::default()
        };

        let indexed = scan_vcf(path.to_str().unwrap(), region.clone()).unwrap();
        assert!(has_tabix_index(path.to_str().unwrap()));
        assert_eq!(indexed.len(), 1);
        assert_eq!(indexed.pos[0], 100);
    }
}