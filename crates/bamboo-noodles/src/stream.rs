use crate::error::NoodlesError;
use crate::lazy_fetch::LazyIndexedFetch;
use crate::record::AlignedRecord;
use bamboo_core::BamScanOptions;
use bamboo_io::{BamSource, BamStorage};
use noodles::bam as bam;
use noodles::bgzf as bgzf;
use noodles::sam::Header;

use std::fs::File;
use std::io::{BufReader, Cursor};
use std::num::NonZeroUsize;
use std::path::Path;

/// A lazy record stream over a BAM source.
pub struct BamRecordStream {
    inner: StreamInner,
}

enum StreamInner {
    LocalScan(LocalScan),
    RemoteScan(RemoteScan),
    LocalFetch(LocalFetch),
    RemoteFetch(RemoteFetch),
}

struct LocalScan {
    reader: bam::io::Reader<bgzf::Reader<File>>,
    header: Header,
    options: BamScanOptions,
    record: bam::Record,
}

struct RemoteScan {
    reader: bam::io::Reader<bgzf::Reader<Cursor<std::sync::Arc<[u8]>>>>,
    header: Header,
    options: BamScanOptions,
    record: bam::Record,
}

struct LocalFetch {
    fetch: LazyIndexedFetch<bgzf::Reader<File>>,
    options: BamScanOptions,
}

struct RemoteFetch {
    fetch: LazyIndexedFetch<bgzf::Reader<Cursor<std::sync::Arc<[u8]>>>>,
    options: BamScanOptions,
}

impl BamRecordStream {
    pub fn open(source: &BamSource, header: &Header, options: BamScanOptions) -> Result<Self, NoodlesError> {
        if options.region.is_some() && bamboo_io::has_index(source) {
            Self::open_fetch(source, options)
        } else {
            Self::open_scan(source, header, options)
        }
    }

    fn open_scan(source: &BamSource, header: &Header, options: BamScanOptions) -> Result<Self, NoodlesError> {
        match &source.storage {
            BamStorage::Local(path) => {
                let mut reader = bam::io::reader::Builder::default()
                    .build_from_path(path)
                    .map_err(NoodlesError::from)?;
                reader.read_header().map_err(NoodlesError::from)?;
                Ok(Self {
                    inner: StreamInner::LocalScan(LocalScan {
                        reader,
                        header: header.clone(),
                        options,
                        record: bam::Record::default(),
                    }),
                })
            }
            BamStorage::Remote { data, .. } => {
                let mut reader = bam::io::Reader::new(Cursor::new(data.clone()));
                reader.read_header().map_err(NoodlesError::from)?;
                Ok(Self {
                    inner: StreamInner::RemoteScan(RemoteScan {
                        reader,
                        header: header.clone(),
                        options,
                        record: bam::Record::default(),
                    }),
                })
            }
        }
    }

    fn open_fetch(source: &BamSource, options: BamScanOptions) -> Result<Self, NoodlesError> {
        let region = options
            .region
            .as_ref()
            .ok_or_else(|| NoodlesError::Message("fetch requires a region".to_string()))?;
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
                let fetch = LazyIndexedFetch::open(reader, header, &parsed_region)?;
                Ok(Self {
                    inner: StreamInner::LocalFetch(LocalFetch { fetch, options }),
                })
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
                let fetch = LazyIndexedFetch::open(reader, header, &parsed_region)?;
                Ok(Self {
                    inner: StreamInner::RemoteFetch(RemoteFetch { fetch, options }),
                })
            }
        }
    }
}

impl Iterator for BamRecordStream {
    type Item = Result<AlignedRecord, NoodlesError>;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.inner {
            StreamInner::LocalScan(state) => next_from_bam_reader(
                &mut state.reader,
                &state.header,
                &state.options,
                &mut state.record,
            ),
            StreamInner::RemoteScan(state) => next_from_bam_reader(
                &mut state.reader,
                &state.header,
                &state.options,
                &mut state.record,
            ),
            StreamInner::LocalFetch(state) => next_from_lazy_fetch(
                &mut state.fetch,
                &state.options,
            ),
            StreamInner::RemoteFetch(state) => next_from_lazy_fetch(
                &mut state.fetch,
                &state.options,
            ),
        }
    }
}

fn next_from_bam_reader<R: std::io::Read>(
    reader: &mut bam::io::Reader<bgzf::Reader<R>>,
    header: &Header,
    options: &BamScanOptions,
    record: &mut bam::Record,
) -> Option<Result<AlignedRecord, NoodlesError>> {
    loop {
        let block_size = match reader.read_record(record) {
            Ok(size) => size,
            Err(err) => return Some(Err(err.into())),
        };
        if block_size == 0 {
            return None;
        }

        let aligned = AlignedRecord::from_bam_record(header, record, options);
        if aligned.passes_filters(options) {
            return Some(Ok(aligned));
        }
    }
}

fn next_from_lazy_fetch<R>(
    fetch: &mut LazyIndexedFetch<R>,
    options: &BamScanOptions,
) -> Option<Result<AlignedRecord, NoodlesError>>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    loop {
        let record = match fetch.next_record() {
            Ok(Some(record)) => record,
            Ok(None) => return None,
            Err(err) => return Some(Err(err)),
        };

        let header = fetch.header();
        let aligned = AlignedRecord::from_bam_record(header, &record, options);
        if aligned.passes_filters(options) {
            return Some(Ok(aligned));
        }
    }
}

pub fn count_records(source: &BamSource, _header: &Header) -> Result<usize, NoodlesError> {
    match &source.storage {
        BamStorage::Local(path) => count_local_path(path),
        BamStorage::Remote { data, .. } => count_remote_bytes(data),
    }
}

fn count_local_path(path: &Path) -> Result<usize, NoodlesError> {
    let workers = std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    let file = File::open(path).map_err(NoodlesError::from)?;
    let decoder = bgzf::MultithreadedReader::with_worker_count(workers, BufReader::new(file));
    let mut reader = bam::io::Reader::from(decoder);
    reader.read_header().map_err(NoodlesError::from)?;
    count_reader_records(&mut reader)
}

fn count_remote_bytes(data: &std::sync::Arc<[u8]>) -> Result<usize, NoodlesError> {
    let mut reader = bam::io::Reader::new(Cursor::new(data.clone()));
    reader.read_header().map_err(NoodlesError::from)?;
    count_reader_records(&mut reader)
}

fn count_reader_records<R: std::io::Read>(
    reader: &mut bam::io::Reader<R>,
) -> Result<usize, NoodlesError> {
    let mut total = 0usize;
    let mut record = bam::Record::default();
    loop {
        let block_size = reader.read_record(&mut record).map_err(NoodlesError::from)?;
        if block_size == 0 {
            break;
        }
        total += 1;
    }
    Ok(total)
}