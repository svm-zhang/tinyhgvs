use pyo3::prelude::*;
use pyo3::types::PyModule;

use tinyhgvs::{ParseHgvsError, ParseHgvsErrorKind};

const PY_ERRORS_MODULE: &str = "tinyhgvs.errors";

pub(crate) struct PyErrorFactory<'py> {
    module: Bound<'py, PyModule>,
}

impl<'py> PyErrorFactory<'py> {
    pub(crate) fn import(py: Python<'py>) -> PyResult<Self> {
        Ok(Self {
            module: PyModule::import(py, PY_ERRORS_MODULE)?,
        })
    }

    fn class(&self, name: &str) -> PyResult<Bound<'py, PyAny>> {
        self.module.getattr(name)
    }

    fn error_kind(&self, value: ParseHgvsErrorKind) -> PyResult<Bound<'py, PyAny>> {
        let value = match value {
            ParseHgvsErrorKind::InvalidSyntax => "invalid_syntax",
            ParseHgvsErrorKind::UnsupportedSyntax => "unsupported_syntax",
            ParseHgvsErrorKind::SemanticConstraint => "semantic_constraint",
        };

        self.class("ParseHgvsErrorKind")?.call1((value,))
    }

    pub(crate) fn parse_error(&self, value: &ParseHgvsError) -> PyResult<PyErr> {
        let error = self.class("TinyHGVSError")?.call1((
            self.error_kind(value.kind())?,
            value.code(),
            value.message(),
            value.input(),
            value.fragment(),
            value.parser_version(),
        ))?;

        Ok(PyErr::from_value(error))
    }
}
