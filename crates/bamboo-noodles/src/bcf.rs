use crate::error::NoodlesError;
use crate::vcf::{append_record_buf, passes_region_filter_buf};
use bamboo_core::{VcfScanOptions, VcfTable};
use noodles::bcf as bcf;
use noodles::csi::{self as csi, binning_index::index::reference_sequence::index::BinnedIndex};
use noodles::vcf as vcf;
use noodles::vcf::variant::Record as _;
use std::ffi::OsString;
use std::io;
use std::path::{Path, PathBuf};

/// High-level BCF reader backed by noodles.
pub struct BcfReader {
    path: String,
    header: vcf::Header,
}

impl BcfReader {
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
        has_csi_index(&self.path)
    }

    pub fn reference_names(&self) -> Vec<String> {
        self.header
            .contigs()
            .keys()
            .map(|name| name.to_string())
            .collect()
    }

    pub fn count_records(&self) -> Result<usize, NoodlesError> {
        scan_bcf(&self.path, VcfScanOptions::default()).map(|table| table.len())
    }

    pub fn scan(&self, options: VcfScanOptions) -> Result<VcfTable, NoodlesError> {
        scan_bcf(&self.path, options)
    }
}

pub fn scan_bcf(path: &str, options: VcfScanOptions) -> Result<VcfTable, NoodlesError> {
    if options.region.is_some() && has_csi_index(path) {
        scan_bcf_indexed(path, &options)
    } else {
        scan_bcf_linear(path, &options)
    }
}

pub fn has_csi_index(path: &str) -> bool {
    csi_index_path(path).exists()
}

/// Build a CSI sidecar for a coordinate-sorted BCF file.
pub fn index_bcf(path: &Path) -> io::Result<csi::Index> {
    use csi::binning_index::{
        index::{
            header::Builder,
            reference_sequence::bin::Chunk,
        },
        Indexer,
    };
    use std::fs::File;

    let file = File::open(path)?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader.read_header()?;

    let mut indexer = Indexer::<BinnedIndex>::default().set_header(Builder::vcf().build());

    let mut record = bcf::Record::default();
    loop {
        let start = reader.get_ref().virtual_position();
        match reader.read_record(&mut record)? {
            0 => break,
            _ => {
                let end = reader.get_ref().virtual_position();
                let reference_sequence_id = record.reference_sequence_id()?;
                let start_position = record
                    .variant_start()
                    .transpose()?
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing position"))?;
                let end_position = record.variant_end(&header)?;
                indexer.add_record(
                    Some((reference_sequence_id, start_position, end_position, true)),
                    Chunk::new(start, end),
                )?;
            }
        }
    }

    Ok(indexer.build(header.contigs().len()))
}

/// Write a CSI sidecar for a BCF file.
pub fn write_bcf_index(path: &Path) -> io::Result<()> {
    let index = index_bcf(path)?;
    let mut index_path = path.as_os_str().to_os_string();
    index_path.push(".csi");
    csi::fs::write(PathBuf::from(index_path), &index)
}

fn csi_index_path(path: &str) -> PathBuf {
    let mut candidate = OsString::from(path);
    candidate.push(".csi");
    PathBuf::from(candidate)
}

fn scan_bcf_linear(path: &str, options: &VcfScanOptions) -> Result<VcfTable, NoodlesError> {
    let mut reader = bcf::io::reader::Builder::default()
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
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

fn scan_bcf_indexed(path: &str, options: &VcfScanOptions) -> Result<VcfTable, NoodlesError> {
    let region = options
        .region
        .as_ref()
        .ok_or_else(|| NoodlesError::Message("indexed BCF scan requires a region".to_string()))?;
    let parsed_region: noodles::core::Region = region
        .to_samtools_region()
        .parse()
        .map_err(|err: noodles::core::region::ParseError| NoodlesError::Message(err.to_string()))?;

    let mut reader = bcf::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(NoodlesError::from)?;
    let header = reader.read_header().map_err(NoodlesError::from)?;
    let mut table = VcfTable::new(options.columns.clone());

    for result in reader
        .query(&header, &parsed_region)
        .map_err(NoodlesError::from)?
    {
        let record = result.map_err(NoodlesError::from)?;
        let record_buf =
            vcf::variant::RecordBuf::try_from_variant_record(&header, &record)
                .map_err(NoodlesError::from)?;
        append_record_buf(&mut table, &record_buf, options);
    }

    Ok(table)
}

fn read_header_from_path(path: &str) -> Result<vcf::Header, NoodlesError> {
    let mut reader = bcf::io::reader::Builder::default()
        .build_from_path(Path::new(path))
        .map_err(NoodlesError::from)?;
    reader.read_header().map_err(NoodlesError::from)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{tiny_bcf_path, tiny_vcf_path, write_tiny_bcf, write_tiny_vcf};
    use bamboo_core::FetchRegion;
    use tempfile::tempdir;

    #[test]
    fn reads_tiny_bcf() {
        let dir = tempdir().unwrap();
        let path = tiny_bcf_path(dir.path());
        write_tiny_bcf(&path).unwrap();

        let reader = BcfReader::open(path.to_str().unwrap()).unwrap();
        assert_eq!(reader.count_records().unwrap(), 2);
        assert_eq!(reader.reference_names(), vec!["chr1", "chr2"]);
    }

    #[test]
    fn bcf_matches_vcf_columns() {
        let dir = tempdir().unwrap();
        let bcf_path = tiny_bcf_path(dir.path());
        let vcf_path = tiny_vcf_path(dir.path());
        write_tiny_bcf(&bcf_path).unwrap();
        write_tiny_vcf(&vcf_path).unwrap();

        let bcf_table = scan_bcf(bcf_path.to_str().unwrap(), VcfScanOptions::default()).unwrap();
        let vcf_table = crate::vcf::scan_vcf(vcf_path.to_str().unwrap(), VcfScanOptions::default()).unwrap();

        assert_eq!(bcf_table.chrom, vcf_table.chrom);
        assert_eq!(bcf_table.pos, vcf_table.pos);
        assert_eq!(bcf_table.reference, vcf_table.reference);
        assert_eq!(bcf_table.alt, vcf_table.alt);
    }

    #[test]
    fn filters_bcf_by_region() {
        let dir = tempdir().unwrap();
        let path = tiny_bcf_path(dir.path());
        write_tiny_bcf(&path).unwrap();

        let options = VcfScanOptions {
            region: Some(FetchRegion {
                reference_name: "chr1".to_string(),
                start: Some(99),
                end: Some(200),
            }),
            ..Default::default()
        };
        let table = scan_bcf(path.to_str().unwrap(), options).unwrap();
        assert_eq!(table.len(), 1);
        assert_eq!(table.chrom[0], "chr1");
        assert_eq!(table.pos[0], 100);
    }

    #[test]
    fn indexed_bcf_region_matches_linear_scan() {
        let dir = tempdir().unwrap();
        let path = tiny_bcf_path(dir.path());
        write_tiny_bcf(&path).unwrap();
        super::write_bcf_index(&path).unwrap();

        let options = VcfScanOptions {
            region: Some(FetchRegion {
                reference_name: "chr1".to_string(),
                start: Some(99),
                end: Some(200),
            }),
            ..Default::default()
        };
        let indexed = scan_bcf(path.to_str().unwrap(), options.clone()).unwrap();
        let linear = scan_bcf_linear(path.to_str().unwrap(), &options).unwrap();
        assert_eq!(indexed.chrom, linear.chrom);
        assert_eq!(indexed.pos, linear.pos);
        assert_eq!(indexed.reference, linear.reference);
        assert_eq!(indexed.alt, linear.alt);
    }
}