use bamboo_core::{BamColumn, BamTable, TagValue, bam_column_name};
use pyo3::prelude::*;
use pyo3::types::PyList;

pub(crate) fn table_to_pyarrow(py: Python<'_>, table: &BamTable) -> PyResult<PyObject> {
    let pa = py.import_bound("pyarrow")?;
    let fields = PyList::empty_bound(py);
    let columns = PyList::empty_bound(py);

    for column in &table.columns {
        let name = bam_column_name(*column);
        let array = column_to_array(py, &pa, *column, table)?;
        let field = pa
            .getattr("field")?
            .call1((name, array.getattr("type")?))?;
        fields.append(field)?;
        columns.append(array)?;
    }

    for tag in &table.tags {
        let array = tag_values_to_array(py, &pa, &tag.values)?;
        let field = pa
            .getattr("field")?
            .call1((tag.name.as_str(), array.getattr("type")?))?;
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

fn column_to_array<'py>(
    py: Python<'py>,
    pa: &Bound<'py, PyModule>,
    column: BamColumn,
    table: &BamTable,
) -> PyResult<Bound<'py, PyAny>> {
    match column {
        BamColumn::QueryName => string_array(py, pa, &table.qname),
        BamColumn::Flag => int_array(py, pa, table.flag.iter().map(|v| Some(*v as i64)).collect()),
        BamColumn::ReferenceName => string_array(py, pa, &table.rname),
        BamColumn::Position => int_array(py, pa, table.pos.iter().map(|v| v.map(i64::from)).collect()),
        BamColumn::MappingQuality => int_array(py, pa, table.mapq.iter().map(|v| v.map(i64::from)).collect()),
        BamColumn::Cigar => {
            let values: Vec<Option<String>> = table.cigar.iter().cloned().map(Some).collect();
            string_array(py, pa, &values)
        }
        BamColumn::MateReferenceName => string_array(py, pa, &table.rnext),
        BamColumn::MatePosition => int_array(py, pa, table.pnext.iter().map(|v| v.map(i64::from)).collect()),
        BamColumn::TemplateLength => int_array(py, pa, table.tlen.iter().map(|v| v.map(i64::from)).collect()),
        BamColumn::Sequence => string_array(py, pa, &table.seq),
        BamColumn::Quality => string_array(py, pa, &table.qual),
    }
}

fn string_array<'py>(
    py: Python<'py>,
    pa: &Bound<'py, PyModule>,
    values: &[Option<String>],
) -> PyResult<Bound<'py, PyAny>> {
    let _ = py;
    pa.getattr("array")?.call1((values.to_vec(),))
}

fn int_array<'py>(
    py: Python<'py>,
    pa: &Bound<'py, PyModule>,
    values: Vec<Option<i64>>,
) -> PyResult<Bound<'py, PyAny>> {
    let _ = py;
    pa.getattr("array")?.call1((values,))
}

fn tag_values_to_array<'py>(
    py: Python<'py>,
    pa: &Bound<'py, PyModule>,
    values: &[TagValue],
) -> PyResult<Bound<'py, PyAny>> {
    let py_values: Vec<PyObject> = values
        .iter()
        .map(|value| match value {
            TagValue::Int(v) => v.to_object(py),
            TagValue::Float(v) => v.to_object(py),
            TagValue::String(v) => v.to_object(py),
            TagValue::Missing => py.None(),
        })
        .collect();
    pa.getattr("array")?.call1((py_values,))
}