use crate::errors;
use crate::{build_scan_options, table::table_to_pyarrow};
use bamboo_core::{BamScanOptions, FetchRegion};
use bamboo_noodles::{scan_bam, AlignedRecord, BamReader};
use pyo3::exceptions::PyStopIteration;
use pyo3::prelude::*;
use pyo3::types::PyDict;


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

#[pyclass(name = "AlignmentIterator")]
pub struct PyAlignmentIterator {
    records: Vec<AlignedRecord>,
    index: usize,
}

#[pymethods]
impl PyAlignmentIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<PyAlignedSegment>> {
        if slf.index >= slf.records.len() {
            return Err(PyStopIteration::new_err(()));
        }
        let record = slf.records[slf.index].clone();
        slf.index += 1;
        Ok(Some(PyAlignedSegment { inner: record }))
    }
}

#[pyclass(name = "AlignmentFile")]
pub struct PyAlignmentFile {
    reader: BamReader,
}

#[pymethods]
impl PyAlignmentFile {
    #[new]
    #[pyo3(signature = (path, mode="rb"))]
    fn new(path: String, mode: &str) -> PyResult<Self> {
        if mode != "rb" && mode != "r" {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Bamboo MVP only supports read mode, got '{mode}'"
            )));
        }
        let reader = BamReader::open(&path).map_err(errors::noodles_to_py_err)?;
        Ok(Self { reader })
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __exit__(
        &self,
        _exc_type: Option<PyObject>,
        _exc: Option<PyObject>,
        _traceback: Option<PyObject>,
    ) -> PyResult<bool> {
        Ok(false)
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<PyAlignmentIterator> {
        slf.fetch(None, None, None, None, None)
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

        let options = BamScanOptions {
            region: fetch_region,
            min_mapq,
            ..BamScanOptions::default()
        };

        let records = self
            .reader
            .iter_records(&options)
            .map_err(errors::noodles_to_py_err)?;
        Ok(PyAlignmentIterator { records, index: 0 })
    }

    fn count(&self) -> PyResult<usize> {
        self.reader.count_records().map_err(errors::noodles_to_py_err)
    }

    fn references(&self) -> Vec<String> {
        self.reader.reference_names()
    }

    #[getter]
    fn reference_lengths(&self) -> Vec<u32> {
        self.reader.reference_lengths()
    }

    fn has_index(&self) -> bool {
        self.reader.has_index()
    }

    fn header(&self, py: Python<'_>) -> PyResult<PyObject> {
        let dict = PyDict::new_bound(py);
        let names = self.reader.reference_names();
        let lengths = self.reader.reference_lengths();
        for (name, length) in names.iter().zip(lengths.iter()) {
            dict.set_item(name, length)?;
        }
        Ok(dict.into())
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
        let table = scan_bam(self.reader.uri(), options)
        .map_err(errors::noodles_to_py_err)?;
        table_to_pyarrow(py, &table)
    }

    fn filename(&self) -> String {
        self.reader.uri().to_string()
    }
}

fn tag_value_to_py(py: Python<'_>, value: &bamboo_core::TagValue) -> PyResult<PyObject> {
    Ok(match value {
        bamboo_core::TagValue::Int(v) => v.to_object(py),
        bamboo_core::TagValue::Float(v) => v.to_object(py),
        bamboo_core::TagValue::String(v) => v.to_object(py),
        bamboo_core::TagValue::Missing => py.None(),
    })
}