mod alignment;
mod cram_file;
mod errors;
mod pileup;
mod table;
mod variant;
mod vcf_table;

use alignment::{PyAlignedSegment, PyAlignmentFile, PyAlignmentIterator};
use cram_file::PyCramFile;
use variant::{PyVariantFile, PyVariantRecord};
use bamboo_core::{BamColumn, BamScanOptions, FetchRegion};
use bamboo_noodles::{scan_bam, scan_cram};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use table::table_to_pyarrow;

#[pymodule]
fn _bamboo(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyCramFile>()?;
    m.add_class::<PyAlignmentFile>()?;
    m.add_class::<PyAlignedSegment>()?;
    m.add_class::<PyAlignmentIterator>()?;
    m.add_function(wrap_pyfunction!(read_bam_table, m)?)?;
    m.add_function(wrap_pyfunction!(scan_bam_table, m)?)?;
    m.add_function(wrap_pyfunction!(read_columns, m)?)?;
    m.add_function(wrap_pyfunction!(read_cram_columns, m)?)?;
    m.add_class::<PyVariantFile>()?;
    m.add_class::<PyVariantRecord>()?;
    m.add_function(wrap_pyfunction!(read_vcf_table, m)?)?;
    m.add_function(wrap_pyfunction!(read_bcf_table, m)?)?;
    m.add_class::<pileup::PyPileupColumn>()?;
    m.add_class::<pileup::PyPileupRead>()?;
    m.add_class::<pileup::PyPileupIterator>()?;
    m.add_function(wrap_pyfunction!(pileup_available, m)?)?;
    m.add_function(wrap_pyfunction!(htslib_available, m)?)?;
    m.add_function(wrap_pyfunction!(primary_backend, m)?)?;
    Ok(())
}

#[pyfunction]
fn pileup_available() -> bool {
    bamboo_htslib::pileup_available()
}

#[pyfunction]
fn htslib_available() -> bool {
    bamboo_htslib::is_available()
}

#[pyfunction]
fn primary_backend() -> &'static str {
    bamboo_htslib::primary_backend()
}

#[pyfunction]
#[pyo3(signature = (path, *, columns=None, region=None))]
fn read_vcf_table(
    py: Python<'_>,
    path: String,
    columns: Option<Vec<String>>,
    region: Option<String>,
) -> PyResult<PyObject> {
    variant::read_variant_table_impl(py, path, columns, region)
}

#[pyfunction]
#[pyo3(signature = (path, *, columns=None, region=None))]
fn read_bcf_table(
    py: Python<'_>,
    path: String,
    columns: Option<Vec<String>>,
    region: Option<String>,
) -> PyResult<PyObject> {
    variant::read_variant_table_impl(py, path, columns, region)
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

#[pyfunction]
#[pyo3(signature = (path, *, columns=None, tags=None, region=None, min_mapq=None, reference_filename=None))]
fn read_cram_columns(
    py: Python<'_>,
    path: String,
    columns: Option<Vec<String>>,
    tags: Option<Vec<String>>,
    region: Option<String>,
    min_mapq: Option<u8>,
    reference_filename: Option<String>,
) -> PyResult<PyObject> {
    let options = build_scan_options(columns, tags, region, min_mapq, None)?;
    let table = scan_cram(&path, options, reference_filename.as_deref())
        .map_err(errors::noodles_to_py_err)?;
    table_to_pyarrow(py, &table)
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