use crate::errors;
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
    fn new(path: String) -> PyResult<Self> {
        let reader = CramReader::open(&path).map_err(errors::noodles_to_py_err)?;
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

    fn count(&self) -> PyResult<usize> {
        self.reader
            .count_records()
            .map_err(errors::noodles_to_py_err)
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
}