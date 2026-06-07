use crate::errors;
use crate::pileup;
use crate::{build_scan_options, table::table_to_pyarrow};
use bamboo_core::{BamScanOptions, FetchRegion};
use bamboo_core::BamTable;
use bamboo_noodles::{
    scan_reader, AlignedRecord, BamReader, BamRecordStream, BamWriter, CramRecordStream,
};
use pyo3::exceptions::{PyStopIteration, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

enum AlignmentFileInner {
    Read(BamReader),
    Write(BamWriter),
}

#[pyclass(name = "AlignedSegment")]
pub struct PyAlignedSegment {
    inner: AlignedRecord,
}

#[pymethods]
impl PyAlignedSegment {
    #[getter]
    fn query_name(&self) -> Option<String> {
        self.inner.query_name.clone()
    }

    #[getter]
    fn flag(&self) -> u16 {
        self.inner.flag
    }

    #[getter]
    fn reference_name(&self) -> Option<String> {
        self.inner.reference_name.clone()
    }

    #[getter]
    fn reference_start(&self) -> Option<i64> {
        self.inner.reference_start
    }

    #[getter]
    fn mapping_quality(&self) -> Option<u8> {
        self.inner.mapping_quality
    }

    #[getter]
    fn cigarstring(&self) -> String {
        self.inner.cigar.clone()
    }

    #[getter]
    fn query_sequence(&self) -> Option<String> {
        self.inner.query_sequence.clone()
    }

    #[getter]
    fn query_qualities(&self) -> Option<String> {
        self.inner.query_qualities.clone()
    }

    fn get_tags(&self, py: Python<'_>) -> PyResult<PyObject> {
        let dict = PyDict::new_bound(py);
        for (name, value) in &self.inner.tags {
            dict.set_item(name, tag_value_to_py(py, value)?)?;
        }
        Ok(dict.into())
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignedSegment(query_name={:?}, reference_name={:?}, reference_start={:?})",
            self.inner.query_name, self.inner.reference_name, self.inner.reference_start
        )
    }
}

enum AlignmentIterInner {
    Stream(BamRecordStream),
    Table { table: BamTable, index: usize },
    CramStream(CramRecordStream),
}

#[pyclass(name = "AlignmentIterator")]
pub struct PyAlignmentIterator {
    inner: AlignmentIterInner,
}

impl PyAlignmentIterator {
    fn from_stream(stream: BamRecordStream) -> Self {
        Self {
            inner: AlignmentIterInner::Stream(stream),
        }
    }

    fn from_table(table: BamTable) -> Self {
        Self {
            inner: AlignmentIterInner::Table { table, index: 0 },
        }
    }

    pub(crate) fn from_cram_stream(stream: CramRecordStream) -> Self {
        Self {
            inner: AlignmentIterInner::CramStream(stream),
        }
    }
}

#[pymethods]
impl PyAlignmentIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<PyAlignedSegment>> {
        match &mut slf.inner {
            AlignmentIterInner::Stream(stream) => match stream.next() {
                Some(Ok(record)) => Ok(Some(PyAlignedSegment { inner: record })),
                Some(Err(err)) => Err(errors::noodles_to_py_err(err)),
                None => Err(PyStopIteration::new_err(())),
            },
            AlignmentIterInner::Table { table, index } => {
                if *index >= table.len() {
                    return Err(PyStopIteration::new_err(()));
                }
                let record = AlignedRecord::from_table_row(table, *index);
                *index += 1;
                Ok(Some(PyAlignedSegment { inner: record }))
            }
            AlignmentIterInner::CramStream(stream) => match stream.next() {
                Some(Ok(record)) => Ok(Some(PyAlignedSegment { inner: record })),
                Some(Err(err)) => Err(errors::noodles_to_py_err(err)),
                None => Err(PyStopIteration::new_err(())),
            },
        }
    }
}

#[pyclass(name = "AlignmentFile")]
pub struct PyAlignmentFile {
    inner: AlignmentFileInner,
}

#[pymethods]
impl PyAlignmentFile {
    #[new]
    #[pyo3(signature = (path, mode="rb", *, header=None, template=None))]
    fn new(
        path: String,
        mode: &str,
        header: Option<&Bound<'_, PyDict>>,
        template: Option<PyRef<'_, PyAlignmentFile>>,
    ) -> PyResult<Self> {
        match mode {
            "rb" | "r" => {
                if header.is_some() || template.is_some() {
                    return Err(PyValueError::new_err(
                        "header and template are only valid for write mode",
                    ));
                }
                let reader = BamReader::open(&path).map_err(errors::noodles_to_py_err)?;
                Ok(Self {
                    inner: AlignmentFileInner::Read(reader),
                })
            }
            "wb" | "w" => {
                let writer = open_write_handle(&path, header, template)?;
                Ok(Self {
                    inner: AlignmentFileInner::Write(writer),
                })
            }
            other => Err(PyValueError::new_err(format!(
                "unsupported AlignmentFile mode '{other}' (use 'rb', 'r', 'wb', or 'w')"
            ))),
        }
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<PyObject>,
        _exc: Option<PyObject>,
        _traceback: Option<PyObject>,
    ) -> PyResult<bool> {
        if let AlignmentFileInner::Write(writer) = &mut self.inner {
            writer.finish().map_err(errors::noodles_to_py_err)?;
        }
        Ok(false)
    }

    fn close(&mut self) -> PyResult<()> {
        if let AlignmentFileInner::Write(writer) = &mut self.inner {
            writer.finish().map_err(errors::noodles_to_py_err)?;
        }
        Ok(())
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<PyAlignmentIterator> {
        let reader = slf.reader()?;
        let options = BamScanOptions::iteration_defaults();
        let stream = reader
            .open_stream(options)
            .map_err(errors::noodles_to_py_err)?;
        Ok(PyAlignmentIterator::from_stream(stream))
    }

    #[pyo3(signature = (contig=None, start=None, stop=None, region=None, min_mapq=None))]
    fn fetch(
        &self,
        contig: Option<String>,
        start: Option<u32>,
        stop: Option<u32>,
        region: Option<String>,
        min_mapq: Option<u8>,
    ) -> PyResult<PyAlignmentIterator> {
        let reader = self.reader().map_err(errors::into_py_err)?;

        let fetch_region = if let Some(region) = region {
            Some(FetchRegion::from_samtools_region(&region).map_err(errors::into_py_err)?)
        } else if let Some(contig) = contig {
            Some(FetchRegion {
                reference_name: contig,
                start,
                end: stop,
            })
        } else {
            None
        };

        Ok(fetch_alignment_iter(reader, fetch_region, min_mapq)?)
    }

    #[pyo3(signature = (contig=None, start=None, stop=None, region=None, min_mapq=None))]
    fn fetch_bulk(
        &self,
        contig: Option<String>,
        start: Option<u32>,
        stop: Option<u32>,
        region: Option<String>,
        min_mapq: Option<u8>,
    ) -> PyResult<Vec<PyAlignedSegment>> {
        let reader = self.reader().map_err(errors::into_py_err)?;

        let fetch_region = if let Some(region) = region {
            Some(FetchRegion::from_samtools_region(&region).map_err(errors::into_py_err)?)
        } else if let Some(contig) = contig {
            Some(FetchRegion {
                reference_name: contig,
                start,
                end: stop,
            })
        } else {
            None
        };

        if fetch_region.is_some() {
            let table = fetch_columnar_table(reader, fetch_region, min_mapq)?;
            return Ok((0..table.len())
                .map(|row| PyAlignedSegment {
                    inner: AlignedRecord::from_table_row(&table, row),
                })
                .collect());
        }

        let mut options = BamScanOptions::iteration_defaults();
        options.min_mapq = min_mapq;

        reader
            .open_stream(options)
            .map_err(errors::noodles_to_py_err)?
            .map(|result| {
                result
                    .map(|record| PyAlignedSegment { inner: record })
                    .map_err(errors::noodles_to_py_err)
            })
            .collect()
    }

    #[pyo3(signature = (*, columns=None, tags=None, region=None, min_mapq=None, reference_name=None))]
    fn fetch_arrow(
        &self,
        py: Python<'_>,
        columns: Option<Vec<String>>,
        tags: Option<Vec<String>>,
        region: Option<String>,
        min_mapq: Option<u8>,
        reference_name: Option<String>,
    ) -> PyResult<PyObject> {
        self.to_arrow(py, columns, tags, region, min_mapq, reference_name)
    }

    fn write(&mut self, record: &PyAlignedSegment) -> PyResult<()> {
        match &mut self.inner {
            AlignmentFileInner::Write(writer) => writer
                .write_record(&record.inner)
                .map_err(errors::noodles_to_py_err),
            AlignmentFileInner::Read(_) => Err(PyValueError::new_err(
                "AlignmentFile is not open for writing",
            )),
        }
    }

    fn count(&self) -> PyResult<usize> {
        self.reader()?
            .count_records()
            .map_err(errors::noodles_to_py_err)
    }

    fn references(&self) -> PyResult<Vec<String>> {
        Ok(self.reader()?.reference_names())
    }

    #[getter]
    fn reference_lengths(&self) -> PyResult<Vec<u32>> {
        Ok(self.reader()?.reference_lengths())
    }

    fn has_index(&self) -> PyResult<bool> {
        Ok(self.reader()?.has_index())
    }

    fn header(&self, py: Python<'_>) -> PyResult<PyObject> {
        let (names, lengths) = match &self.inner {
            AlignmentFileInner::Read(reader) => {
                (reader.reference_names(), reader.reference_lengths())
            }
            AlignmentFileInner::Write(writer) => {
                let header = writer.header();
                let names: Vec<String> = header
                    .reference_sequences()
                    .keys()
                    .map(|name| name.to_string())
                    .collect();
                let lengths: Vec<u32> = header
                    .reference_sequences()
                    .values()
                    .map(|reference| reference.length().get() as u32)
                    .collect();
                (names, lengths)
            }
        };

        let dict = PyDict::new_bound(py);
        for (name, length) in names.iter().zip(lengths.iter()) {
            dict.set_item(name, length)?;
        }
        Ok(dict.into())
    }

    #[pyo3(signature = (*, columns=None, tags=None, region=None, min_mapq=None, reference_name=None))]
    fn read_table(
        &self,
        py: Python<'_>,
        columns: Option<Vec<String>>,
        tags: Option<Vec<String>>,
        region: Option<String>,
        min_mapq: Option<u8>,
        reference_name: Option<String>,
    ) -> PyResult<PyObject> {
        self.to_arrow(py, columns, tags, region, min_mapq, reference_name)
    }

    #[pyo3(signature = (*, columns=None, tags=None, region=None, min_mapq=None, reference_name=None))]
    fn to_arrow(
        &self,
        py: Python<'_>,
        columns: Option<Vec<String>>,
        tags: Option<Vec<String>>,
        region: Option<String>,
        min_mapq: Option<u8>,
        reference_name: Option<String>,
    ) -> PyResult<PyObject> {
        let options = build_scan_options(columns, tags, region, min_mapq, reference_name)?;
        let table = match &self.inner {
            AlignmentFileInner::Read(reader) => {
                scan_reader(reader, options).map_err(errors::noodles_to_py_err)?
            }
            AlignmentFileInner::Write(_) => {
                return Err(PyValueError::new_err(
                    "to_arrow() requires a read-mode AlignmentFile",
                ));
            }
        };
        table_to_pyarrow(py, &table)
    }

    #[pyo3(signature = (contig=None, start=0, stop=0x7fff_ffff, *, region=None, reads=true))]
    fn pileup(
        &self,
        contig: Option<String>,
        start: u32,
        stop: u32,
        region: Option<String>,
        reads: bool,
    ) -> PyResult<pileup::PyPileupIterator> {
        let reader = self.reader()?;
        if !reader.has_index() {
            return Err(PyValueError::new_err(
                "pileup requires an indexed BAM (missing .bai)",
            ));
        }

        let (reference_name, region_start, region_end) = if let Some(region) = region {
            let parsed =
                FetchRegion::from_samtools_region(&region).map_err(errors::into_py_err)?;
            let start = parsed.start.unwrap_or(0);
            let end = parsed.end.unwrap_or(u32::MAX);
            (parsed.reference_name, start, end)
        } else {
            let contig = contig.ok_or_else(|| {
                PyValueError::new_err("pileup requires contig or region=...")
            })?;
            (contig, start, stop)
        };

        pileup::pileup_region(
            reader.uri(),
            &reference_name,
            region_start,
            region_end,
            None,
            reads,
        )
    }

    fn filename(&self) -> String {
        match &self.inner {
            AlignmentFileInner::Read(reader) => reader.uri().to_string(),
            AlignmentFileInner::Write(writer) => writer.path().to_string(),
        }
    }

    #[getter]
    fn mode(&self) -> &'static str {
        match self.inner {
            AlignmentFileInner::Read(_) => "rb",
            AlignmentFileInner::Write(_) => "wb",
        }
    }
}

impl PyAlignmentFile {
    fn reader(&self) -> Result<&BamReader, PyErr> {
        match &self.inner {
            AlignmentFileInner::Read(reader) => Ok(reader),
            AlignmentFileInner::Write(_) => Err(PyValueError::new_err(
                "this operation requires a read-mode AlignmentFile",
            )),
        }
    }
}

fn open_write_handle(
    path: &str,
    header: Option<&Bound<'_, PyDict>>,
    template: Option<PyRef<'_, PyAlignmentFile>>,
) -> PyResult<BamWriter> {
    if let Some(template) = template {
        return match &template.inner {
            AlignmentFileInner::Read(reader) => {
                BamWriter::create_from_reader(path, reader).map_err(errors::noodles_to_py_err)
            }
            AlignmentFileInner::Write(_) => Err(PyValueError::new_err(
                "template must be a read-mode AlignmentFile",
            )),
        };
    }

    let header_dict = header.ok_or_else(|| {
        PyValueError::new_err("write mode requires header= or template=")
    })?;

    let mut refs = Vec::new();
    for (key, value) in header_dict.iter() {
        let name: String = key.extract()?;
        let length: u32 = value.extract()?;
        refs.push((name, length));
    }

    BamWriter::create_with_references(path, &refs).map_err(errors::noodles_to_py_err)
}

fn fetch_alignment_iter(
    reader: &BamReader,
    fetch_region: Option<FetchRegion>,
    min_mapq: Option<u8>,
) -> PyResult<PyAlignmentIterator> {
    if fetch_region.is_some() {
        let table = fetch_columnar_table(reader, fetch_region, min_mapq)?;
        return Ok(PyAlignmentIterator::from_table(table));
    }

    let mut options = BamScanOptions::iteration_defaults();
    options.min_mapq = min_mapq;
    let stream = reader
        .open_stream(options)
        .map_err(errors::noodles_to_py_err)?;
    Ok(PyAlignmentIterator::from_stream(stream))
}

fn fetch_columnar_table(
    reader: &BamReader,
    fetch_region: Option<FetchRegion>,
    min_mapq: Option<u8>,
) -> PyResult<BamTable> {
    let mut options = BamScanOptions::iteration_defaults();
    options.region = fetch_region;
    options.min_mapq = min_mapq;
    scan_reader(reader, options).map_err(errors::noodles_to_py_err)
}

fn tag_value_to_py(py: Python<'_>, value: &bamboo_core::TagValue) -> PyResult<PyObject> {
    Ok(match value {
        bamboo_core::TagValue::Int(v) => v.to_object(py),
        bamboo_core::TagValue::Float(v) => v.to_object(py),
        bamboo_core::TagValue::String(v) => v.to_object(py),
        bamboo_core::TagValue::Missing => py.None(),
    })
}