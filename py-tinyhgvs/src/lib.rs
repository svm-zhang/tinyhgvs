mod error;
mod model_codec;

use pyo3::prelude::*;
use pyo3::types::PyModule;
use tinyhgvs::parse_hgvs as parse_hgvs_core;

use crate::error::PyErrorFactory;
use crate::model_codec::PyModelCodec;

#[pyfunction]
fn parse_hgvs<'py>(py: Python<'py>, input: &str) -> PyResult<Bound<'py, PyAny>> {
    match parse_hgvs_core(input) {
        Ok(variant) => PyModelCodec::import(py)?.hgvs_variant(&variant),
        Err(error) => Err(PyErrorFactory::import(py)?.parse_error(&error)?),
    }
}

#[pymodule]
fn _tinyhgvs(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_hgvs, m)?)?;
    Ok(())
}
