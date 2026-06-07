use pyo3::exceptions::{PyRuntimeError, PyStopIteration, PyValueError};
use pyo3::prelude::*;

#[pyclass(name = "PileupRead")]
#[derive(Clone)]
pub struct PyPileupRead {
    pub query_name: Option<String>,
    pub query_position: Option<u32>,
    pub is_del: bool,
    pub is_head: bool,
    pub is_tail: bool,
    pub is_refskip: bool,
}

#[pymethods]
impl PyPileupRead {
    #[getter]
    fn query_name(&self) -> Option<String> {
        self.query_name.clone()
    }

    #[getter]
    fn query_position(&self) -> Option<u32> {
        self.query_position
    }

    #[getter]
    fn is_del(&self) -> bool {
        self.is_del
    }

    #[getter]
    fn is_head(&self) -> bool {
        self.is_head
    }

    #[getter]
    fn is_tail(&self) -> bool {
        self.is_tail
    }

    #[getter]
    fn is_refskip(&self) -> bool {
        self.is_refskip
    }

    fn __repr__(&self) -> String {
        format!(
            "PileupRead(query_name={:?}, query_position={:?})",
            self.query_name, self.query_position
        )
    }
}

#[pyclass(name = "PileupColumn")]
#[derive(Clone)]
pub struct PyPileupColumn {
    pub reference_name: String,
    pub reference_id: i32,
    pub position: u32,
    pub depth: u32,
    pub reads: Vec<PyPileupRead>,
}

#[pymethods]
impl PyPileupColumn {
    #[getter]
    fn reference_name(&self) -> &str {
        &self.reference_name
    }

    #[getter]
    fn reference_id(&self) -> i32 {
        self.reference_id
    }

    #[getter]
    fn pos(&self) -> u32 {
        self.position
    }

    #[getter]
    fn n(&self) -> u32 {
        self.depth
    }

    #[getter]
    fn pileups(&self) -> Vec<PyPileupRead> {
        self.reads.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "PileupColumn(reference_name={:?}, pos={}, n={})",
            self.reference_name, self.position, self.depth
        )
    }
}

#[pyclass(name = "PileupIterator", unsendable)]
pub struct PyPileupIterator {
    #[cfg(feature = "htslib")]
    stream: Option<bamboo_htslib::PileupStream>,
}

#[pymethods]
impl PyPileupIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<PyPileupColumn>> {
        #[cfg(feature = "htslib")]
        {
            let stream = slf
                .stream
                .as_mut()
                .ok_or_else(|| PyRuntimeError::new_err("pileup iterator is exhausted"))?;

            match stream.next_column() {
                Some(Ok(column)) => Ok(Some(column_to_py(&column))),
                Some(Err(err)) => Err(htslib_to_py_err(err)),
                None => {
                    slf.stream = None;
                    Err(PyStopIteration::new_err(()))
                }
            }
        }

        #[cfg(not(feature = "htslib"))]
        {
            Err(PyRuntimeError::new_err(
                "pileup requires Bamboo built with the 'htslib' feature",
            ))
        }
    }
}

pub fn pileup_region(
    path: &str,
    contig: &str,
    start: u32,
    end: u32,
    reference_filename: Option<&str>,
    materialize_reads: bool,
) -> PyResult<PyPileupIterator> {
    #[cfg(feature = "htslib")]
    {
        let stream = bamboo_htslib::PileupStream::open(
            path,
            contig,
            start,
            end,
            reference_filename,
            materialize_reads,
        )
        .map_err(htslib_to_py_err)?;

        return Ok(PyPileupIterator {
            stream: Some(stream),
        });
    }

    #[cfg(not(feature = "htslib"))]
    {
        let _ = (
            path,
            contig,
            start,
            end,
            reference_filename,
            materialize_reads,
        );
        Err(PyRuntimeError::new_err(
            "pileup requires Bamboo built with the 'htslib' feature",
        ))
    }
}

#[cfg(feature = "htslib")]
fn column_to_py(column: &bamboo_htslib::PileupColumn) -> PyPileupColumn {
    PyPileupColumn {
        reference_name: column.reference_name.clone(),
        reference_id: column.reference_id,
        position: column.position,
        depth: column.depth,
        reads: column
            .reads
            .iter()
            .map(|read| PyPileupRead {
                query_name: read.query_name.clone(),
                query_position: read.query_position,
                is_del: read.is_del,
                is_head: read.is_head,
                is_tail: read.is_tail,
                is_refskip: read.is_refskip,
            })
            .collect(),
    }
}

#[cfg(feature = "htslib")]
fn htslib_to_py_err(error: bamboo_htslib::HtslibError) -> PyErr {
    match error {
        bamboo_htslib::HtslibError::Io(err) => PyValueError::new_err(err.to_string()),
        bamboo_htslib::HtslibError::Htslib(err) => PyRuntimeError::new_err(err.to_string()),
        bamboo_htslib::HtslibError::Message(message) => PyRuntimeError::new_err(message),
    }
}