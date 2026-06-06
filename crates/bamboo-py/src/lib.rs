mod alignment;
mod errors;
mod table;

use alignment::{PyAlignedSegment, PyAlignmentFile, PyAlignmentIterator};
use bamboo_core::{BamColumn, BamScanOptions, FetchRegion};
use bamboo_noodles::scan_bam;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use table::table_to_pyarrow;

#[pymodule]
fn _bamboo(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyAlignmentFile>()?;
    m.add_class::<PyAlignedSegment>()?;
    m.add_class::<PyAlignmentIterator>()?;
    m.add_function(wrap_pyfunction!(read_bam_table, m)?)?;
    m.add_function(wrap_pyfunction!(scan_bam_table, m)?)?;
    m.add_function(wrap_pyfunction!(read_columns, m)?)?;
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (path, *, columns=None, tags=None, region=None, min_mapq=None, reference_name=None))]
fn read_bam_table(
    py: Python<'_>,
    path: String,
    columns: Option<Vec<String>>,
    tags: Option<Vec<String>>,
    region: Option<String>,
    min_mapq: Option<u8>,
    reference_name: Option<String>,
) -> PyResult<PyObject> {
    let options = build_scan_options(columns, tags, region, min_mapq, reference_name)?;
    let table = scan_bam(&path, options).map_err(errors::noodles_to_py_err)?;
    table_to_pyarrow(py, &table)
}

#[pyfunction]
#[pyo3(signature = (path, *, columns=None, tags=None, region=None, min_mapq=None, reference_name=None))]
fn scan_bam_table(
    py: Python<'_>,
    path: String,
    columns: Option<Vec<String>>,
    tags: Option<Vec<String>>,
    region: Option<String>,
    min_mapq: Option<u8>,
    reference_name: Option<String>,
) -> PyResult<PyObject> {
    read_bam_table(
        py,
        path,
        columns,
        tags,
        region,
        min_mapq,
        reference_name,
    )
}

/// Fast columnar scan — the recommended production path for analytics workloads.
#[pyfunction]
#[pyo3(signature = (path, *, columns=None, tags=None, region=None, min_mapq=None, reference_name=None))]
fn read_columns(
    py: Python<'_>,
    path: String,
    columns: Option<Vec<String>>,
    tags: Option<Vec<String>>,
    region: Option<String>,
    min_mapq: Option<u8>,
    reference_name: Option<String>,
) -> PyResult<PyObject> {
    read_bam_table(
        py,
        path,
        columns,
        tags,
        region,
        min_mapq,
        reference_name,
    )
}

pub(crate) fn build_scan_options(
    columns: Option<Vec<String>>,
    tags: Option<Vec<String>>,
    region: Option<String>,
    min_mapq: Option<u8>,
    reference_name: Option<String>,
) -> PyResult<BamScanOptions> {
    let columns = match columns {
        Some(names) => names
            .iter()
            .map(|name| {
                BamColumn::parse_name(name).ok_or_else(|| {
                    PyValueError::new_err(format!("unknown BAM column '{name}'"))
                })
            })
            .collect::<PyResult<Vec<_>>>()?,
        None => BamScanOptions::default().columns,
    };

    let region = match region {
        Some(value) => Some(FetchRegion::from_samtools_region(&value).map_err(errors::into_py_err)?),
        None => None,
    };

    Ok(BamScanOptions {
        columns,
        tags: tags.unwrap_or_default(),
        region,
        min_mapq,
        reference_name,
    })
}