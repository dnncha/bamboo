use crate::alignment::PyAlignmentIterator;
use crate::errors;
use crate::pileup;
use crate::{build_scan_options, table::table_to_pyarrow};
use bamboo_core::{BamScanOptions, FetchRegion};
use bamboo_noodles::CramReader;
use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass(name = "CramFile")]
pub struct PyCramFile {
    reader: CramReader,
}

#[pymethods]
impl PyCramFile {
    #[new]
    #[pyo3(signature = (path, *, reference_filename=None))]
    fn new(path: String, reference_filename: Option<String>) -> PyResult<Self> {
        let reader = CramReader::open_with_reference(
            &path,
            reference_filename.as_deref(),
        )
        .map_err(errors::noodles_to_py_err)?;
        Ok(Self { reader })
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
        Ok(false)
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<PyAlignmentIterator> {
        let options = BamScanOptions::iteration_defaults();
        let stream = slf
            .reader
            .open_stream(options)
            .map_err(errors::noodles_to_py_err)?;
        Ok(PyAlignmentIterator::from_cram_stream(stream))
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

        let mut options = BamScanOptions::iteration_defaults();
        options.region = fetch_region;
        options.min_mapq = min_mapq;

        let stream = self
            .reader
            .open_stream(options)
            .map_err(errors::noodles_to_py_err)?;
        Ok(PyAlignmentIterator::from_cram_stream(stream))
    }

    fn count(&self) -> PyResult<usize> {
        self.reader
            .count_records()
            .map_err(errors::noodles_to_py_err)
    }

    fn has_index(&self) -> bool {
        self.reader.has_index()
    }

    fn references(&self) -> Vec<String> {
        self.reader.reference_names()
    }

    #[getter]
    fn reference_lengths(&self) -> Vec<u32> {
        self.reader.reference_lengths()
    }

    fn filename(&self) -> String {
        self.reader.path().to_string()
    }

    fn header(&self, py: Python<'_>) -> PyResult<PyObject> {
        let names = self.reader.reference_names();
        let lengths = self.reader.reference_lengths();
        let dict = PyDict::new_bound(py);
        for (name, length) in names.iter().zip(lengths.iter()) {
            dict.set_item(name, length)?;
        }
        Ok(dict.into())
    }

    #[pyo3(signature = (*, columns=None, tags=None, region=None, min_mapq=None))]
    fn read_table(
        &self,
        py: Python<'_>,
        columns: Option<Vec<String>>,
        tags: Option<Vec<String>>,
        region: Option<String>,
        min_mapq: Option<u8>,
    ) -> PyResult<PyObject> {
        self.to_arrow(py, columns, tags, region, min_mapq)
    }

    #[pyo3(signature = (contig=None, start=0, stop=0x7fff_ffff, *, region=None))]
    fn pileup(
        &self,
        contig: Option<String>,
        start: u32,
        stop: u32,
        region: Option<String>,
    ) -> PyResult<pileup::PyPileupIterator> {
        if !self.reader.has_index() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "pileup requires an indexed CRAM (missing .crai)",
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
                pyo3::exceptions::PyValueError::new_err("pileup requires contig or region=...")
            })?;
            (contig, start, stop)
        };

        pileup::pileup_region(
            self.reader.path(),
            &reference_name,
            region_start,
            region_end,
            self.reader.reference_fasta(),
        )
    }

    #[pyo3(signature = (*, columns=None, tags=None, region=None, min_mapq=None))]
    fn to_arrow(
        &self,
        py: Python<'_>,
        columns: Option<Vec<String>>,
        tags: Option<Vec<String>>,
        region: Option<String>,
        min_mapq: Option<u8>,
    ) -> PyResult<PyObject> {
        let options = build_scan_options(columns, tags, region, min_mapq, None)?;
        let table = self
            .reader
            .scan(options)
            .map_err(errors::noodles_to_py_err)?;
        table_to_pyarrow(py, &table)
    }
}