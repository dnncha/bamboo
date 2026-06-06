use pyo3::prelude::*;

/// Bamboo - High-performance HTS data for Python
///
/// This is the Rust core (exposed via PyO3) for the Bamboo library.
/// The goal is to eventually provide fast, safe, Arrow-native access to
/// BAM/CRAM, VCF/BCF and related formats, with excellent cloud + streaming support.

#[pymodule]
fn _bamboo(py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("hello", py_fn!(py, hello()))?;

    // Future:
    // m.add_class::<AlignmentFile>()?;
    // m.add_class::<VariantFile>()?;
    // etc.

    Ok(())
}

#[pyfunction]
fn hello() -> String {
    "Hello from Bamboo (Rust + PyO3)! The modern pysam successor is waking up.".to_string()
}

// Placeholder for future high-level API
// We will build nice Rust types here that are then exposed cleanly to Python,
// with heavy use of noodles for format parsing and arrow-rs for the DataFrame story.
