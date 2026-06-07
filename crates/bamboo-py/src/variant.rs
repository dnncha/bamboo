use crate::errors;
use crate::vcf_table::vcf_table_to_pyarrow;
use bamboo_core::{FetchRegion, VcfColumn, VcfScanOptions};
use bamboo_noodles::{scan_bcf, scan_vcf, VariantReader};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass(name = "VariantRecord")]
pub struct PyVariantRecord {
    chrom: String,
    pos: i32,
    id: String,
    reference: String,
    alt: String,
    qual: Option<f32>,
    filter: String,
}

#[pymethods]
impl PyVariantRecord {
    #[getter]
    fn chrom(&self) -> &str {
        &self.chrom
    }

    #[getter]
    fn pos(&self) -> i32 {
        self.pos
    }

    #[getter]
    fn id(&self) -> &str {
        &self.id
    }

    #[getter]
    fn ref_(&self) -> &str {
        &self.reference
    }

    #[getter]
    fn alt(&self) -> &str {
        &self.alt
    }

    #[getter]
    fn qual(&self) -> Option<f32> {
        self.qual
    }

    #[getter]
    fn filter(&self) -> &str {
        &self.filter
    }

    fn __repr__(&self) -> String {
        format!(
            "VariantRecord(chrom={:?}, pos={}, ref={:?}, alt={:?})",
            self.chrom, self.pos, self.reference, self.alt
        )
    }
}

#[pyclass(name = "VariantFile")]
pub struct PyVariantFile {
    reader: VariantReader,
}

#[pymethods]
impl PyVariantFile {
    #[new]
    fn new(path: String) -> PyResult<Self> {
        let reader = VariantReader::open(&path).map_err(errors::noodles_to_py_err)?;
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
        self.reader.count_records().map_err(errors::noodles_to_py_err)
    }

    fn references(&self) -> Vec<String> {
        self.reader.reference_names()
    }

    fn has_index(&self) -> bool {
        self.reader.has_index()
    }

    fn filename(&self) -> String {
        self.reader.path().to_string()
    }

    fn header(&self, py: Python<'_>) -> PyResult<PyObject> {
        let dict = PyDict::new_bound(py);
        for name in self.reader.reference_names() {
            dict.set_item(name, py.None())?;
        }
        Ok(dict.into())
    }

    #[pyo3(signature = (*, columns=None, region=None))]
    fn to_arrow(
        &self,
        py: Python<'_>,
        columns: Option<Vec<String>>,
        region: Option<String>,
    ) -> PyResult<PyObject> {
        let options = build_vcf_scan_options(columns, region)?;
        let table = self
            .reader
            .scan(options)
            .map_err(errors::noodles_to_py_err)?;
        vcf_table_to_pyarrow(py, &table)
    }

    #[pyo3(signature = (contig=None, start=None, stop=None, region=None))]
    fn fetch(
        &self,
        contig: Option<String>,
        start: Option<u32>,
        stop: Option<u32>,
        region: Option<String>,
    ) -> PyResult<Vec<PyVariantRecord>> {
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

        let options = VcfScanOptions {
            region: fetch_region,
            ..Default::default()
        };
        let table = self
            .reader
            .scan(options)
            .map_err(errors::noodles_to_py_err)?;

        let mut records = Vec::with_capacity(table.len());
        for index in 0..table.len() {
            records.push(PyVariantRecord {
                chrom: table.chrom[index].clone(),
                pos: table.pos[index],
                id: table.id[index].clone(),
                reference: table.reference[index].clone(),
                alt: table.alt[index].clone(),
                qual: table.qual[index],
                filter: table.filter[index].clone(),
            });
        }
        Ok(records)
    }
}

pub(crate) fn read_variant_table_impl(
    py: Python<'_>,
    path: String,
    columns: Option<Vec<String>>,
    region: Option<String>,
) -> PyResult<PyObject> {
    let options = build_vcf_scan_options(columns, region)?;
    let table = if path
        .rsplit_once('.')
        .is_some_and(|(_, extension)| extension.eq_ignore_ascii_case("bcf"))
    {
        scan_bcf(&path, options).map_err(errors::noodles_to_py_err)?
    } else {
        scan_vcf(&path, options).map_err(errors::noodles_to_py_err)?
    };
    vcf_table_to_pyarrow(py, &table)
}

pub(crate) fn build_vcf_scan_options(
    columns: Option<Vec<String>>,
    region: Option<String>,
) -> PyResult<VcfScanOptions> {
    let columns = match columns {
        Some(names) => names
            .iter()
            .map(|name| {
                VcfColumn::parse_name(name).ok_or_else(|| {
                    PyValueError::new_err(format!("unknown VCF column '{name}'"))
                })
            })
            .collect::<PyResult<Vec<_>>>()?,
        None => VcfScanOptions::default().columns,
    };

    let region = match region {
        Some(value) => Some(FetchRegion::from_samtools_region(&value).map_err(errors::into_py_err)?),
        None => None,
    };

    Ok(VcfScanOptions { columns, region })
}