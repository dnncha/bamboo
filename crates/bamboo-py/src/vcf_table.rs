use bamboo_core::{VcfColumn, VcfTable};
use pyo3::prelude::*;
use pyo3::types::PyList;

pub(crate) fn vcf_table_to_pyarrow(py: Python<'_>, table: &VcfTable) -> PyResult<PyObject> {
    let pa = py.import_bound("pyarrow")?;
    let fields = PyList::empty_bound(py);
    let columns = PyList::empty_bound(py);

    for column in &table.columns {
        let name = column.arrow_name();
        let array = vcf_column_to_array(py, &pa, *column, table)?;
        let field = pa
            .getattr("field")?
            .call1((name, array.getattr("type")?))?;
        fields.append(field)?;
        columns.append(array)?;
    }

    let schema = pa.getattr("schema")?.call1((fields,))?;
    let kwargs = pyo3::types::PyDict::new_bound(py);
    kwargs.set_item("schema", schema)?;
    pa.getattr("Table")?
        .call_method("from_arrays", (columns,), Some(&kwargs))
        .map(|value| value.into())
}

fn vcf_column_to_array<'py>(
    py: Python<'py>,
    pa: &Bound<'py, PyModule>,
    column: VcfColumn,
    table: &VcfTable,
) -> PyResult<Bound<'py, PyAny>> {
    match column {
        VcfColumn::Chrom => string_array(py, pa, &table.chrom),
        VcfColumn::Pos => int_array(py, pa, table.pos.iter().map(|v| Some(i64::from(*v))).collect()),
        VcfColumn::Id => string_array(py, pa, &table.id),
        VcfColumn::Ref => string_array(py, pa, &table.reference),
        VcfColumn::Alt => string_array(py, pa, &table.alt),
        VcfColumn::Qual => float_array(py, pa, table.qual.iter().map(|v| v.map(f64::from)).collect()),
        VcfColumn::Filter => string_array(py, pa, &table.filter),
    }
}

fn string_array<'py>(
    py: Python<'py>,
    pa: &Bound<'py, PyModule>,
    values: &[String],
) -> PyResult<Bound<'py, PyAny>> {
    let _ = py;
    let optional: Vec<Option<String>> = values.iter().cloned().map(Some).collect();
    pa.getattr("array")?.call1((optional,))
}

fn int_array<'py>(
    py: Python<'py>,
    pa: &Bound<'py, PyModule>,
    values: Vec<Option<i64>>,
) -> PyResult<Bound<'py, PyAny>> {
    let _ = py;
    pa.getattr("array")?.call1((values,))
}

fn float_array<'py>(
    py: Python<'py>,
    pa: &Bound<'py, PyModule>,
    values: Vec<Option<f64>>,
) -> PyResult<Bound<'py, PyAny>> {
    let _ = py;
    pa.getattr("array")?.call1((values,))
}