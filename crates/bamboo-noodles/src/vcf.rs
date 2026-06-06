use crate::error::NoodlesError;
use bamboo_core::{FetchRegion, VcfColumn, VcfScanOptions, VcfTable};
use noodles::vcf as vcf;
use noodles::vcf::variant::record::Ids;
use std::io::BufRead;
use std::path::Path;

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
    let mut reader = open_reader(path)?;
    let header = reader.read_header().map_err(NoodlesError::from)?;
    let mut table = VcfTable::new(options.columns.clone());

    for result in reader.record_bufs(&header) {
        let record = result.map_err(NoodlesError::from)?;
        if !passes_region_filter(&record, &options.region) {
            continue;
        }
        append_record(&mut table, &record, &options);
    }

    Ok(table)
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

fn passes_region_filter(record: &vcf::variant::RecordBuf, region: &Option<FetchRegion>) -> bool {
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

fn append_record(table: &mut VcfTable, record: &vcf::variant::RecordBuf, options: &VcfScanOptions) {
    let chrom = record.reference_sequence_name().to_string();
    let pos = record
        .variant_start()
        .map(|position| position.get() as i32)
        .unwrap_or(0);
    let id = record
        .ids()
        .iter()
        .next()
        .unwrap_or(".")
        .to_string();
    let reference = record.reference_bases().to_string();
    let alt = record
        .alternate_bases()
        .as_ref()
        .iter()
        .map(|allele| allele.as_str())
        .collect::<Vec<_>>()
        .join(",");
    let qual = record.quality_score();
    let filter = record
        .filters()
        .as_ref()
        .iter()
        .next()
        .cloned()
        .unwrap_or_else(|| ".".to_string());

    if options.wants_column(VcfColumn::Chrom) {
        table.chrom.push(chrom);
    }
    if options.wants_column(VcfColumn::Pos) {
        table.pos.push(pos);
    }
    if options.wants_column(VcfColumn::Id) {
        table.id.push(id);
    }
    if options.wants_column(VcfColumn::Ref) {
        table.reference.push(reference);
    }
    if options.wants_column(VcfColumn::Alt) {
        table.alt.push(alt);
    }
    if options.wants_column(VcfColumn::Qual) {
        table.qual.push(qual);
    }
    if options.wants_column(VcfColumn::Filter) {
        table.filter.push(filter);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{tiny_vcf_path, write_tiny_vcf};
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
}