use pyo3::exceptions::{PyRuntimeError, PyStopIteration, PyValueError};
use pyo3::ffi;
use pyo3::prelude::*;
use pyo3::types::PyString;
use std::cell::Cell;
use std::sync::{Arc, OnceLock};

#[cfg(feature = "htslib")]
fn next_column_without_gil(
    stream: &mut bamboo_htslib::PileupStream,
) -> Option<Result<bamboo_htslib::PileupColumn, bamboo_htslib::HtslibError>> {
    // htslib pileup iterators hold raw C pointers (not Send); release the GIL via the C-API.
    unsafe {
        let thread_state = ffi::PyEval_SaveThread();
        let result = stream.next_column();
        ffi::PyEval_RestoreThread(thread_state);
        result
    }
}

fn empty_reads() -> Arc<Vec<PyPileupRead>> {
    static EMPTY: OnceLock<Arc<Vec<PyPileupRead>>> = OnceLock::new();
    Arc::clone(EMPTY.get_or_init(|| Arc::new(Vec::new())))
}

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

/// Thin pileup column proxy modeled after pysam's PileupColumn.
///
/// When the pileup iterator runs with `reads=False`, a single column instance is
/// updated in place each step (same model as pysam's live pileup buffer).
#[pyclass(name = "PileupColumn", freelist = 2048)]
#[derive(Clone)]
pub struct PyPileupColumn {
    reference_names: Arc<Vec<Py<PyString>>>,
    fallback_contig: Arc<Py<PyString>>,
    reference_id: Cell<i32>,
    position: Cell<u32>,
    depth: Cell<u32>,
    reads: Arc<Vec<PyPileupRead>>,
}

#[pymethods]
impl PyPileupColumn {
    #[getter]
    fn reference_name(&self, py: Python<'_>) -> Py<PyString> {
        self.reference_names
            .get(self.reference_id.get() as usize)
            .cloned()
            .unwrap_or_else(|| self.fallback_contig.clone_ref(py).clone())
    }

    #[getter]
    fn reference_id(&self) -> i32 {
        self.reference_id.get()
    }

    #[getter]
    fn pos(&self) -> u32 {
        self.position.get()
    }

    #[getter]
    fn n(&self) -> u32 {
        self.depth.get()
    }

    #[getter]
    fn pileups(&self) -> Vec<PyPileupRead> {
        self.reads.as_ref().clone()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "PileupColumn(reference_name={:?}, pos={}, n={})",
            self.reference_name(py).to_string(),
            self.position.get(),
            self.depth.get()
        ))
    }
}

#[pyclass(name = "PileupIterator", unsendable)]
pub struct PyPileupIterator {
    #[cfg(feature = "htslib")]
    stream: Option<bamboo_htslib::PileupStream>,
    reference_names: Arc<Vec<Py<PyString>>>,
    fallback_contig: Arc<Py<PyString>>,
    materialize_reads: bool,
    /// Reused and mutated in place when `reads=False` to avoid per-column allocation.
    count_proxy: Option<Py<PyPileupColumn>>,
}

#[pymethods]
impl PyPileupIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Py<PyPileupColumn>> {
        #[cfg(feature = "htslib")]
        {
            let py = slf.py();
            let stream = slf
                .stream
                .as_mut()
                .ok_or_else(|| PyRuntimeError::new_err("pileup iterator is exhausted"))?;

            let next = next_column_without_gil(stream);
            match next {
                Some(Ok(column)) => {
                    if slf.materialize_reads {
                        let reads = Arc::new(
                            column
                                .reads
                                .into_iter()
                                .map(|read| PyPileupRead {
                                    query_name: read.query_name,
                                    query_position: read.query_position,
                                    is_del: read.is_del,
                                    is_head: read.is_head,
                                    is_tail: read.is_tail,
                                    is_refskip: read.is_refskip,
                                })
                                .collect(),
                        );
                        return Ok(Py::new(
                            py,
                            PyPileupColumn {
                                reference_names: Arc::clone(&slf.reference_names),
                                fallback_contig: Arc::clone(&slf.fallback_contig),
                                reference_id: Cell::new(column.reference_id),
                                position: Cell::new(column.position),
                                depth: Cell::new(column.depth),
                                reads,
                            },
                        )?);
                    }

                    if slf.count_proxy.is_none() {
                        slf.count_proxy = Some(Py::new(
                            py,
                            PyPileupColumn {
                                reference_names: Arc::clone(&slf.reference_names),
                                fallback_contig: Arc::clone(&slf.fallback_contig),
                                reference_id: Cell::new(column.reference_id),
                                position: Cell::new(column.position),
                                depth: Cell::new(column.depth),
                                reads: empty_reads(),
                            },
                        )?);
                    }

                    let proxy = slf.count_proxy.as_ref().unwrap();
                    let current = proxy.borrow(py);
                    current.reference_id.set(column.reference_id);
                    current.position.set(column.position);
                    current.depth.set(column.depth);
                    Ok(proxy.clone_ref(py))
                }
                Some(Err(err)) => Err(htslib_to_py_err(err)),
                None => {
                    slf.stream = None;
                    slf.count_proxy = None;
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
    py: Python<'_>,
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

        let reference_names = Arc::new(
            stream
                .target_names()
                .iter()
                .map(|name| PyString::new_bound(py, name).unbind())
                .collect::<Vec<_>>(),
        );
        let fallback_contig = Arc::new(PyString::new_bound(py, contig).unbind());

        return Ok(PyPileupIterator {
            stream: Some(stream),
            reference_names,
            fallback_contig,
            materialize_reads,
            count_proxy: None,
        });
    }

    #[cfg(not(feature = "htslib"))]
    {
        let _ = (
            py,
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
fn htslib_to_py_err(error: bamboo_htslib::HtslibError) -> PyErr {
    match error {
        bamboo_htslib::HtslibError::Io(err) => PyValueError::new_err(err.to_string()),
        bamboo_htslib::HtslibError::Htslib(err) => PyRuntimeError::new_err(err.to_string()),
        bamboo_htslib::HtslibError::Message(message) => PyRuntimeError::new_err(message),
    }
}