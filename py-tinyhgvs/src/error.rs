const PY_ERRORS_MODULE: &str = "tinyhgvs.errors";
const PY_MODELS_MODULE: &str = "tinyhgvs.models";

struct PyErrorFactory<'py> {
    module: Bound<'py, PyModule>,
}

impl<'py> PyErrorFactory<'py> {
    fn import(py: Python<'py>) -> PyResult<Self> {
        Ok(Self {
            module: PyModule::import(py, PY_ERRORS_MODULE)?,
        })
    }

    fn class(&self, name: &str) -> PyResult<Bound<'py, PyAny>> {
        self.module.getattr(name)
    }

    fn error_kind(&self, value: ParseHgvsErrorKind) -> PyResult<Bound<'py, PyAny>> {
        self.class("ParseHgvsErrorKind")?
            .call1((error_kind_value(value),))
    }

    fn parse_error(&self, value: &CoreParseHgvsError) -> PyResult<PyErr> {
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

fn error_kind_value(value: ParseHgvsErrorKind) -> &'static str {
    match value {
        ParseHgvsErrorKind::InvalidSyntax => "invalid_syntax",
        ParseHgvsErrorKind::UnsupportedSyntax => "unsupported_syntax",
        ParseHgvsErrorKind::SemanticConstraint => "semantic_constraint",
    }
}
