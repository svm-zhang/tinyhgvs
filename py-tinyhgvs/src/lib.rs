use pyo3::prelude::*;
use pyo3::types::{PyModule, PyTuple};
use tinyhgvs::{
    parse_hgvs as parse_hgvs_core, Accession, CoordinateSystem, CopiedSequenceItem,
    HgvsVariant as CoreHgvsVariant, Interval, LiteralSequenceItem, NucleotideAnchor,
    NucleotideCoordinate, NucleotideEdit, NucleotideSequenceItem,
    ParseHgvsError as CoreParseHgvsError, ParseHgvsErrorKind, ProteinCoordinate, ProteinEdit,
    ProteinEffect, ProteinSequence, ProteinVariant, ReferenceSpec, RepeatSequenceItem,
    VariantDescription,
};

const PY_ERRORS_MODULE: &str = "tinyhgvs.errors";
const PY_MODELS_MODULE: &str = "tinyhgvs.models";

#[pyfunction]
fn parse_hgvs<'py>(py: Python<'py>, input: &str) -> PyResult<Bound<'py, PyAny>> {
    match parse_hgvs_core(input) {
        Ok(variant) => PyModelCodec::import(py)?.variant(&variant),
        Err(error) => Err(PyErrorFactory::import(py)?.parse_error(&error)?),
    }
}

#[pyfunction]
fn _roundtrip_variant<'py>(
    py: Python<'py>,
    value: Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    let codec = PyModelCodec::import(py)?;
    let variant = codec.extract_variant(&value)?;
    codec.variant(&variant)
}

#[pymodule]
fn _tinyhgvs(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_hgvs, m)?)?;
    m.add_function(wrap_pyfunction!(_roundtrip_variant, m)?)?;
    Ok(())
}

struct PyModelCodec<'py> {
    py: Python<'py>,
    module: Bound<'py, PyModule>,
}

impl<'py> PyModelCodec<'py> {
    fn import(py: Python<'py>) -> PyResult<Self> {
        Ok(Self {
            py,
            module: PyModule::import(py, PY_MODELS_MODULE)?,
        })
    }

    fn class(&self, name: &str) -> PyResult<Bound<'py, PyAny>> {
        self.module.getattr(name)
    }

    fn class_name(&self, value: &Bound<'py, PyAny>) -> PyResult<String> {
        value.getattr("__class__")?.getattr("__name__")?.extract()
    }

    fn enum_value(&self, value: &Bound<'py, PyAny>) -> PyResult<String> {
        match value.getattr("value") {
            Ok(raw) => raw.extract(),
            Err(_) => value.extract(),
        }
    }

    fn coordinate_system(&self, value: CoordinateSystem) -> PyResult<Bound<'py, PyAny>> {
        self.class("CoordinateSystem")?.call1((value.as_str(),))
    }

    fn coordinate_system_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<CoordinateSystem> {
        match self.enum_value(value)?.as_str() {
            "g" => Ok(CoordinateSystem::Genomic),
            "o" => Ok(CoordinateSystem::CircularGenomic),
            "m" => Ok(CoordinateSystem::Mitochondrial),
            "c" => Ok(CoordinateSystem::CodingDna),
            "n" => Ok(CoordinateSystem::NonCodingDna),
            "r" => Ok(CoordinateSystem::Rna),
            "p" => Ok(CoordinateSystem::Protein),
            other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "unsupported CoordinateSystem value: {other}"
            ))),
        }
    }

    fn nucleotide_anchor(&self, value: NucleotideAnchor) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotideAnchor")?
            .call1((nucleotide_anchor_value(value),))
    }

    fn nucleotide_anchor_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<NucleotideAnchor> {
        match self.enum_value(value)?.as_str() {
            "absolute" => Ok(NucleotideAnchor::Absolute),
            "relative_cds_start" => Ok(NucleotideAnchor::RelativeCdsStart),
            "relative_cds_end" => Ok(NucleotideAnchor::RelativeCdsEnd),
            other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "unsupported NucleotideAnchor value: {other}"
            ))),
        }
    }

    fn accession(&self, value: &Accession) -> PyResult<Bound<'py, PyAny>> {
        self.class("Accession")?.call1((&value.id, value.version))
    }

    fn accession_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<Accession> {
        Ok(Accession {
            id: value.getattr("id")?.extract()?,
            version: value.getattr("version")?.extract()?,
        })
    }

    fn reference_spec(&self, value: &ReferenceSpec) -> PyResult<Bound<'py, PyAny>> {
        let context = value
            .context
            .as_ref()
            .map(|context| self.accession(context))
            .transpose()?;

        self.class("ReferenceSpec")?
            .call1((self.accession(&value.primary)?, context))
    }

    fn reference_spec_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<ReferenceSpec> {
        let context = value.getattr("context")?;
        Ok(ReferenceSpec {
            primary: self.accession_from_py(&value.getattr("primary")?)?,
            context: if context.is_none() {
                None
            } else {
                Some(self.accession_from_py(&context)?)
            },
        })
    }

    fn nucleotide_coordinate(&self, value: &NucleotideCoordinate) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotideCoordinate")?.call1((
            self.nucleotide_anchor(value.anchor)?,
            value.coordinate,
            value.offset,
        ))
    }

    fn nucleotide_coordinate_from_py(
        &self,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<NucleotideCoordinate> {
        Ok(NucleotideCoordinate {
            anchor: self.nucleotide_anchor_from_py(&value.getattr("anchor")?)?,
            coordinate: value.getattr("coordinate")?.extract()?,
            offset: value.getattr("offset")?.extract()?,
        })
    }

    fn protein_coordinate(&self, value: &ProteinCoordinate) -> PyResult<Bound<'py, PyAny>> {
        self.class("ProteinCoordinate")?
            .call1((&value.residue, value.ordinal))
    }

    fn protein_coordinate_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<ProteinCoordinate> {
        Ok(ProteinCoordinate {
            residue: value.getattr("residue")?.extract()?,
            ordinal: value.getattr("ordinal")?.extract()?,
        })
    }

    fn nucleotide_interval(
        &self,
        value: &Interval<NucleotideCoordinate>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let end = value
            .end
            .as_ref()
            .map(|end| self.nucleotide_coordinate(end))
            .transpose()?;

        self.class("Interval")?
            .call1((self.nucleotide_coordinate(&value.start)?, end))
    }

    fn nucleotide_interval_from_py(
        &self,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<Interval<NucleotideCoordinate>> {
        let end = value.getattr("end")?;
        Ok(Interval {
            start: self.nucleotide_coordinate_from_py(&value.getattr("start")?)?,
            end: if end.is_none() {
                None
            } else {
                Some(self.nucleotide_coordinate_from_py(&end)?)
            },
        })
    }

    fn protein_interval(&self, value: &Interval<ProteinCoordinate>) -> PyResult<Bound<'py, PyAny>> {
        let end = value
            .end
            .as_ref()
            .map(|end| self.protein_coordinate(end))
            .transpose()?;

        self.class("Interval")?
            .call1((self.protein_coordinate(&value.start)?, end))
    }

    fn protein_interval_from_py(
        &self,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<Interval<ProteinCoordinate>> {
        let end = value.getattr("end")?;
        Ok(Interval {
            start: self.protein_coordinate_from_py(&value.getattr("start")?)?,
            end: if end.is_none() {
                None
            } else {
                Some(self.protein_coordinate_from_py(&end)?)
            },
        })
    }

    fn copied_sequence_item(&self, value: &CopiedSequenceItem) -> PyResult<Bound<'py, PyAny>> {
        let source_reference = value
            .source_reference
            .as_ref()
            .map(|reference| self.reference_spec(reference))
            .transpose()?;
        let source_coordinate_system = value
            .source_coordinate_system
            .map(|coordinate_system| self.coordinate_system(coordinate_system))
            .transpose()?;

        self.class("CopiedSequenceItem")?.call1((
            source_reference,
            source_coordinate_system,
            self.nucleotide_interval(&value.source_location)?,
            value.is_inverted,
        ))
    }

    fn copied_sequence_item_from_py(
        &self,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<CopiedSequenceItem> {
        let source_reference = value.getattr("source_reference")?;
        let source_coordinate_system = value.getattr("source_coordinate_system")?;
        Ok(CopiedSequenceItem {
            source_reference: if source_reference.is_none() {
                None
            } else {
                Some(self.reference_spec_from_py(&source_reference)?)
            },
            source_coordinate_system: if source_coordinate_system.is_none() {
                None
            } else {
                Some(self.coordinate_system_from_py(&source_coordinate_system)?)
            },
            source_location: self
                .nucleotide_interval_from_py(&value.getattr("source_location")?)?,
            is_inverted: value.getattr("is_inverted")?.extract()?,
        })
    }

    fn nucleotide_sequence_item(
        &self,
        value: &NucleotideSequenceItem,
    ) -> PyResult<Bound<'py, PyAny>> {
        match value {
            NucleotideSequenceItem::Literal(LiteralSequenceItem { value }) => {
                self.class("LiteralSequenceItem")?.call1((value,))
            }
            NucleotideSequenceItem::Repeat(RepeatSequenceItem { unit, count }) => {
                self.class("RepeatSequenceItem")?.call1((unit, *count))
            }
            NucleotideSequenceItem::Copied(item) => self.copied_sequence_item(item),
        }
    }

    fn nucleotide_sequence_item_from_py(
        &self,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<NucleotideSequenceItem> {
        match self.class_name(value)?.as_str() {
            "LiteralSequenceItem" => Ok(NucleotideSequenceItem::Literal(LiteralSequenceItem {
                value: value.getattr("value")?.extract()?,
            })),
            "RepeatSequenceItem" => Ok(NucleotideSequenceItem::Repeat(RepeatSequenceItem {
                unit: value.getattr("unit")?.extract()?,
                count: value.getattr("count")?.extract()?,
            })),
            "CopiedSequenceItem" => Ok(NucleotideSequenceItem::Copied(
                self.copied_sequence_item_from_py(value)?,
            )),
            other => Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(format!(
                "unsupported nucleotide sequence item type: {other}"
            ))),
        }
    }

    fn nucleotide_items_tuple(
        &self,
        value: &[NucleotideSequenceItem],
    ) -> PyResult<Bound<'py, PyTuple>> {
        let items = value
            .iter()
            .map(|item| self.nucleotide_sequence_item(item))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(self.py, items)
    }

    fn nucleotide_items_from_py(
        &self,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<Vec<NucleotideSequenceItem>> {
        let mut items = Vec::new();
        for item in value.try_iter()? {
            items.push(self.nucleotide_sequence_item_from_py(&item?)?);
        }
        Ok(items)
    }

    fn protein_sequence(&self, value: &ProteinSequence) -> PyResult<Bound<'py, PyAny>> {
        let residues = PyTuple::new(self.py, &value.residues)?;
        self.class("ProteinSequence")?.call1((residues,))
    }

    fn protein_sequence_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<ProteinSequence> {
        Ok(ProteinSequence {
            residues: value.getattr("residues")?.extract()?,
        })
    }

    fn nucleotide_sequence_omitted_edit(
        &self,
        value: &NucleotideEdit,
    ) -> PyResult<Bound<'py, PyAny>> {
        let name = match value {
            NucleotideEdit::NoChange => "no_change",
            NucleotideEdit::Deletion => "deletion",
            NucleotideEdit::Duplication => "duplication",
            NucleotideEdit::Inversion => "inversion",
            _ => unreachable!("nucleotide_sequence_omitted_edit called for payload edit"),
        };
        self.class("NucleotideSequenceOmittedEdit")?.call1((name,))
    }

    fn nucleotide_sequence_omitted_edit_from_py(
        &self,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<NucleotideEdit> {
        match self.enum_value(value)?.as_str() {
            "no_change" => Ok(NucleotideEdit::NoChange),
            "deletion" => Ok(NucleotideEdit::Deletion),
            "duplication" => Ok(NucleotideEdit::Duplication),
            "inversion" => Ok(NucleotideEdit::Inversion),
            other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "unsupported NucleotideSequenceOmittedEdit value: {other}"
            ))),
        }
    }

    fn nucleotide_edit(&self, value: &NucleotideEdit) -> PyResult<Bound<'py, PyAny>> {
        match value {
            NucleotideEdit::NoChange => self.nucleotide_sequence_omitted_edit(value),
            NucleotideEdit::Substitution {
                reference,
                alternate,
            } => self
                .class("NucleotideSubstitutionEdit")?
                .call1((reference, alternate)),
            NucleotideEdit::Deletion => self.nucleotide_sequence_omitted_edit(value),
            NucleotideEdit::Duplication => self.nucleotide_sequence_omitted_edit(value),
            NucleotideEdit::Insertion { items } => self
                .class("NucleotideInsertionEdit")?
                .call1((self.nucleotide_items_tuple(items)?,)),
            NucleotideEdit::Inversion => self.nucleotide_sequence_omitted_edit(value),
            NucleotideEdit::DeletionInsertion { items } => self
                .class("NucleotideDeletionInsertionEdit")?
                .call1((self.nucleotide_items_tuple(items)?,)),
        }
    }

    fn nucleotide_edit_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<NucleotideEdit> {
        match self.class_name(value)?.as_str() {
            "NucleotideSequenceOmittedEdit" => self.nucleotide_sequence_omitted_edit_from_py(value),
            "NucleotideSubstitutionEdit" => Ok(NucleotideEdit::Substitution {
                reference: value.getattr("reference")?.extract()?,
                alternate: value.getattr("alternate")?.extract()?,
            }),
            "NucleotideInsertionEdit" => Ok(NucleotideEdit::Insertion {
                items: self.nucleotide_items_from_py(&value.getattr("items")?)?,
            }),
            "NucleotideDeletionInsertionEdit" => Ok(NucleotideEdit::DeletionInsertion {
                items: self.nucleotide_items_from_py(&value.getattr("items")?)?,
            }),
            other => Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(format!(
                "unsupported nucleotide edit type: {other}"
            ))),
        }
    }

    fn protein_sequence_omitted_edit(&self, value: &ProteinEdit) -> PyResult<Bound<'py, PyAny>> {
        let name = match value {
            ProteinEdit::Unknown => "unknown",
            ProteinEdit::NoChange => "no_change",
            ProteinEdit::Deletion => "deletion",
            ProteinEdit::Duplication => "duplication",
            _ => unreachable!("protein_sequence_omitted_edit called for payload edit"),
        };
        self.class("ProteinSequenceOmittedEdit")?.call1((name,))
    }

    fn protein_sequence_omitted_edit_from_py(
        &self,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<ProteinEdit> {
        match self.enum_value(value)?.as_str() {
            "unknown" => Ok(ProteinEdit::Unknown),
            "no_change" => Ok(ProteinEdit::NoChange),
            "deletion" => Ok(ProteinEdit::Deletion),
            "duplication" => Ok(ProteinEdit::Duplication),
            other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "unsupported ProteinSequenceOmittedEdit value: {other}"
            ))),
        }
    }

    fn protein_edit(&self, value: &ProteinEdit) -> PyResult<Bound<'py, PyAny>> {
        match value {
            ProteinEdit::Unknown => self.protein_sequence_omitted_edit(value),
            ProteinEdit::NoChange => self.protein_sequence_omitted_edit(value),
            ProteinEdit::Substitution { to } => self.class("ProteinSubstitutionEdit")?.call1((to,)),
            ProteinEdit::Deletion => self.protein_sequence_omitted_edit(value),
            ProteinEdit::Duplication => self.protein_sequence_omitted_edit(value),
            ProteinEdit::Insertion { sequence } => self
                .class("ProteinInsertionEdit")?
                .call1((self.protein_sequence(sequence)?,)),
            ProteinEdit::DeletionInsertion { sequence } => self
                .class("ProteinDeletionInsertionEdit")?
                .call1((self.protein_sequence(sequence)?,)),
        }
    }

    fn protein_edit_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<ProteinEdit> {
        match self.class_name(value)?.as_str() {
            "ProteinSequenceOmittedEdit" => self.protein_sequence_omitted_edit_from_py(value),
            "ProteinSubstitutionEdit" => Ok(ProteinEdit::Substitution {
                to: value.getattr("to")?.extract()?,
            }),
            "ProteinInsertionEdit" => Ok(ProteinEdit::Insertion {
                sequence: self.protein_sequence_from_py(&value.getattr("sequence")?)?,
            }),
            "ProteinDeletionInsertionEdit" => Ok(ProteinEdit::DeletionInsertion {
                sequence: self.protein_sequence_from_py(&value.getattr("sequence")?)?,
            }),
            other => Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(format!(
                "unsupported protein edit type: {other}"
            ))),
        }
    }

    fn protein_effect(&self, value: &ProteinEffect) -> PyResult<Bound<'py, PyAny>> {
        match value {
            ProteinEffect::Unknown => self.class("ProteinUnknownEffect")?.call0(),
            ProteinEffect::NoProteinProduced => {
                self.class("ProteinNoProteinProducedEffect")?.call0()
            }
            ProteinEffect::Edit { location, edit } => self
                .class("ProteinEditEffect")?
                .call1((self.protein_interval(location)?, self.protein_edit(edit)?)),
        }
    }

    fn protein_effect_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<ProteinEffect> {
        match self.class_name(value)?.as_str() {
            "ProteinUnknownEffect" => Ok(ProteinEffect::Unknown),
            "ProteinNoProteinProducedEffect" => Ok(ProteinEffect::NoProteinProduced),
            "ProteinEditEffect" => Ok(ProteinEffect::Edit {
                location: self.protein_interval_from_py(&value.getattr("location")?)?,
                edit: self.protein_edit_from_py(&value.getattr("edit")?)?,
            }),
            other => Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(format!(
                "unsupported protein effect type: {other}"
            ))),
        }
    }

    fn description(&self, value: &VariantDescription) -> PyResult<Bound<'py, PyAny>> {
        match value {
            VariantDescription::Nucleotide(value) => self.class("NucleotideVariant")?.call1((
                self.nucleotide_interval(&value.location)?,
                self.nucleotide_edit(&value.edit)?,
            )),
            VariantDescription::Protein(value) => self
                .class("ProteinVariant")?
                .call1((value.is_predicted, self.protein_effect(&value.effect)?)),
        }
    }

    fn description_from_py(&self, value: &Bound<'py, PyAny>) -> PyResult<VariantDescription> {
        match self.class_name(value)?.as_str() {
            "NucleotideVariant" => Ok(VariantDescription::Nucleotide(
                tinyhgvs::NucleotideVariant {
                    location: self.nucleotide_interval_from_py(&value.getattr("location")?)?,
                    edit: self.nucleotide_edit_from_py(&value.getattr("edit")?)?,
                },
            )),
            "ProteinVariant" => Ok(VariantDescription::Protein(ProteinVariant {
                is_predicted: value.getattr("is_predicted")?.extract()?,
                effect: self.protein_effect_from_py(&value.getattr("effect")?)?,
            })),
            other => Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(format!(
                "unsupported variant description type: {other}"
            ))),
        }
    }

    fn variant(&self, value: &CoreHgvsVariant) -> PyResult<Bound<'py, PyAny>> {
        let reference = value
            .reference
            .as_ref()
            .map(|reference| self.reference_spec(reference))
            .transpose()?;

        self.class("HgvsVariant")?.call1((
            reference,
            self.coordinate_system(value.coordinate_system)?,
            self.description(&value.description)?,
        ))
    }

    fn extract_variant(&self, value: &Bound<'py, PyAny>) -> PyResult<CoreHgvsVariant> {
        let reference = value.getattr("reference")?;
        Ok(CoreHgvsVariant {
            reference: if reference.is_none() {
                None
            } else {
                Some(self.reference_spec_from_py(&reference)?)
            },
            coordinate_system: self
                .coordinate_system_from_py(&value.getattr("coordinate_system")?)?,
            description: self.description_from_py(&value.getattr("description")?)?,
        })
    }
}

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

fn nucleotide_anchor_value(value: NucleotideAnchor) -> &'static str {
    match value {
        NucleotideAnchor::Absolute => "absolute",
        NucleotideAnchor::RelativeCdsStart => "relative_cds_start",
        NucleotideAnchor::RelativeCdsEnd => "relative_cds_end",
    }
}

fn error_kind_value(value: ParseHgvsErrorKind) -> &'static str {
    match value {
        ParseHgvsErrorKind::InvalidSyntax => "invalid_syntax",
        ParseHgvsErrorKind::UnsupportedSyntax => "unsupported_syntax",
        ParseHgvsErrorKind::SemanticConstraint => "semantic_constraint",
    }
}
