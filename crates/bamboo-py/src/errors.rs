use bamboo_core::RegionParseError;
use bamboo_noodles::NoodlesError;
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::PyErr;

pub fn into_py_err<E: std::error::Error>(error: E) -> PyErr {
    PyRuntimeError::new_err(error.to_string())
}

pub fn noodles_to_py_err(error: NoodlesError) -> PyErr {
    match error {
        NoodlesError::Io(err) => PyIOError::new_err(err.to_string()),
        NoodlesError::ObjectStore(err) => PyIOError::new_err(err.to_string()),
        NoodlesError::Region(err) => region_to_py_err(err),
        NoodlesError::MissingIndex { path } => {
            PyIOError::new_err(format!("missing BAM index for indexed fetch: {path}"))
        }
        NoodlesError::MissingReference { name } => {
            PyValueError::new_err(format!("reference sequence not found: {name}"))
        }
        NoodlesError::Message(message) => PyRuntimeError::new_err(message),
    }
}

fn region_to_py_err(error: RegionParseError) -> PyErr {
    PyValueError::new_err(error.to_string())
}