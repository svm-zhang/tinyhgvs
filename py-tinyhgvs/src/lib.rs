use pyo3::prelude::*;
use pyo3::types::{PyModule, PyTuple};
use tinyhgvs::{
    parse_hgvs as parse_hgvs_core, Accession, Allele, AllelePhase, AlleleVariant, CoordinateSystem,
    CopiedSequenceItem, HgvsVariant as CoreHgvsVariant, Interval, LiteralSequenceItem,
    NucleotideAnchor, NucleotideCoordinate, NucleotideEdit, NucleotideRepeatBlock,
    NucleotideSequenceItem, NucleotideVariant, ParseHgvsError as CoreParseHgvsError,
    ParseHgvsErrorKind, ProteinCoordinate, ProteinEdit, ProteinEffect, ProteinExtensionEdit,
    ProteinExtensionTerminal, ProteinFrameshiftStop, ProteinFrameshiftStopKind, ProteinSequence,
    ReferenceSpec, RepeatSequenceItem, VariantDescription,
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

#[pymodule]
fn _tinyhgvs(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_hgvs, m)?)?;
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

    fn coordinate_system(&self, value: CoordinateSystem) -> PyResult<Bound<'py, PyAny>> {
        self.class("CoordinateSystem")?.call1((value.as_str(),))
    }

    fn nucleotide_anchor(&self, value: NucleotideAnchor) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotideAnchor")?
            .call1((nucleotide_anchor_value(value),))
    }

    fn accession(&self, value: &Accession) -> PyResult<Bound<'py, PyAny>> {
        self.class("Accession")?.call1((&value.id, value.version))
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

    fn nucleotide_coordinate(&self, value: &NucleotideCoordinate) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotideCoordinate")?.call1((
            self.nucleotide_anchor(value.anchor)?,
            value.coordinate,
            value.offset,
        ))
    }

    fn protein_coordinate(&self, value: &ProteinCoordinate) -> PyResult<Bound<'py, PyAny>> {
        self.class("ProteinCoordinate")?
            .call1((&value.residue, value.ordinal))
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

    fn protein_interval(&self, value: &Interval<ProteinCoordinate>) -> PyResult<Bound<'py, PyAny>> {
        let end = value
            .end
            .as_ref()
            .map(|end| self.protein_coordinate(end))
            .transpose()?;

        self.class("Interval")?
            .call1((self.protein_coordinate(&value.start)?, end))
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

    fn nucleotide_repeat_block(
        &self,
        value: &NucleotideRepeatBlock,
    ) -> PyResult<Bound<'py, PyAny>> {
        let location = value
            .location
            .as_ref()
            .map(|location| self.nucleotide_interval(location))
            .transpose()?;

        self.class("NucleotideRepeatBlock")?
            .call1((value.count, value.unit.as_deref(), location))
    }

    fn nucleotide_repeat_blocks_tuple(
        &self,
        value: &[NucleotideRepeatBlock],
    ) -> PyResult<Bound<'py, PyTuple>> {
        let items = value
            .iter()
            .map(|item| self.nucleotide_repeat_block(item))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(self.py, items)
    }

    fn allele_phase(&self, value: AllelePhase) -> PyResult<Bound<'py, PyAny>> {
        let name = match value {
            AllelePhase::Trans => "trans",
            AllelePhase::Uncertain => "uncertain",
        };
        self.class("AllelePhase")?.call1((name,))
    }

    fn allele<T>(
        &self,
        value: &Allele<T>,
        map_variant: fn(&Self, &T) -> PyResult<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let variants = value
            .variants
            .iter()
            .map(|item| map_variant(self, item))
            .collect::<PyResult<Vec<_>>>()?;

        self.class("Allele")?
            .call1((PyTuple::new(self.py, variants)?,))
    }

    fn alleles_tuple<T>(
        &self,
        value: &[Allele<T>],
        map_variant: fn(&Self, &T) -> PyResult<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let alleles = value
            .iter()
            .map(|item| self.allele(item, map_variant))
            .collect::<PyResult<Vec<_>>>()?;

        PyTuple::new(self.py, alleles)
    }

    fn allele_variant<T>(
        &self,
        value: &AlleleVariant<T>,
        map_variant: fn(&Self, &T) -> PyResult<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let allele_two = value
            .allele_two
            .as_ref()
            .map(|allele| self.allele(allele, map_variant))
            .transpose()?;
        let phase = value
            .phase
            .map(|phase| self.allele_phase(phase))
            .transpose()?;

        // Python reuses the same Allele / AlleleVariant container types for
        // nucleotide and protein allele descriptions.
        self.class("AlleleVariant")?.call1((
            self.allele(&value.allele_one, map_variant)?,
            allele_two,
            phase,
            self.alleles_tuple(&value.alleles_unphased, map_variant)?,
        ))
    }

    fn nucleotide_variant(&self, value: &NucleotideVariant) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotideVariant")?.call1((
            self.nucleotide_interval(&value.location)?,
            self.nucleotide_edit(&value.edit)?,
        ))
    }

    fn protein_variant(&self, value: &tinyhgvs::ProteinVariant) -> PyResult<Bound<'py, PyAny>> {
        self.class("ProteinVariant")?
            .call1((value.is_predicted, self.protein_effect(&value.effect)?))
    }

    fn protein_sequence(&self, value: &ProteinSequence) -> PyResult<Bound<'py, PyAny>> {
        let residues = PyTuple::new(self.py, &value.residues)?;
        self.class("ProteinSequence")?.call1((residues,))
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
            _ => unreachable!("nucleotide_sequence_omitted_edit called for non-enum edit"),
        };
        self.class("NucleotideSequenceOmittedEdit")?.call1((name,))
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
            NucleotideEdit::Repeat { blocks } => self
                .class("NucleotideRepeatEdit")?
                .call1((self.nucleotide_repeat_blocks_tuple(blocks)?,)),
        }
    }

    fn protein_sequence_omitted_edit(&self, value: &ProteinEdit) -> PyResult<Bound<'py, PyAny>> {
        let name = match value {
            ProteinEdit::Unknown => "unknown",
            ProteinEdit::NoChange => "no_change",
            ProteinEdit::Deletion => "deletion",
            ProteinEdit::Duplication => "duplication",
            _ => unreachable!("protein_sequence_omitted_edit called for non-enum edit"),
        };
        self.class("ProteinSequenceOmittedEdit")?.call1((name,))
    }

    fn protein_frameshift_stop_kind(
        &self,
        value: ProteinFrameshiftStopKind,
    ) -> PyResult<Bound<'py, PyAny>> {
        let name = match value {
            ProteinFrameshiftStopKind::Omitted => "omitted",
            ProteinFrameshiftStopKind::Unknown => "unknown",
            ProteinFrameshiftStopKind::Known => "known",
        };
        self.class("ProteinFrameshiftStopKind")?.call1((name,))
    }

    fn protein_frameshift_stop(
        &self,
        value: &ProteinFrameshiftStop,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.class("ProteinFrameshiftStop")?.call1((
            value.ordinal,
            self.protein_frameshift_stop_kind(value.kind)?,
        ))
    }

    fn protein_extension_terminal(
        &self,
        value: ProteinExtensionTerminal,
    ) -> PyResult<Bound<'py, PyAny>> {
        let name = match value {
            ProteinExtensionTerminal::N => "n",
            ProteinExtensionTerminal::C => "c",
        };
        self.class("ProteinExtensionTerminal")?.call1((name,))
    }

    fn protein_extension_edit(&self, value: &ProteinExtensionEdit) -> PyResult<Bound<'py, PyAny>> {
        self.class("ProteinExtensionEdit")?.call1((
            self.protein_extension_terminal(value.to_terminal)?,
            value.to_residue.as_deref(),
            value.terminal_ordinal,
        ))
    }

    fn protein_edit(&self, value: &ProteinEdit) -> PyResult<Bound<'py, PyAny>> {
        match value {
            ProteinEdit::Unknown => self.protein_sequence_omitted_edit(value),
            ProteinEdit::NoChange => self.protein_sequence_omitted_edit(value),
            ProteinEdit::Substitution { to } => self.class("ProteinSubstitutionEdit")?.call1((to,)),
            ProteinEdit::Deletion => self.protein_sequence_omitted_edit(value),
            ProteinEdit::Duplication => self.protein_sequence_omitted_edit(value),
            ProteinEdit::Extension(extension) => self.protein_extension_edit(extension),
            ProteinEdit::Frameshift { to_residue, stop } => self
                .class("ProteinFrameshiftEdit")?
                .call1((to_residue, self.protein_frameshift_stop(stop)?)),
            ProteinEdit::Insertion { sequence } => self
                .class("ProteinInsertionEdit")?
                .call1((self.protein_sequence(sequence)?,)),
            ProteinEdit::DeletionInsertion { sequence } => self
                .class("ProteinDeletionInsertionEdit")?
                .call1((self.protein_sequence(sequence)?,)),
            ProteinEdit::Repeat { count } => self.class("ProteinRepeatEdit")?.call1((count,)),
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

    fn description(&self, value: &VariantDescription) -> PyResult<Bound<'py, PyAny>> {
        match value {
            VariantDescription::Nucleotide(value) => self.nucleotide_variant(value),
            VariantDescription::NucleotideAllele(value) => {
                self.allele_variant(value, Self::nucleotide_variant)
            }
            VariantDescription::Protein(value) => self.protein_variant(value),
            VariantDescription::ProteinAllele(value) => {
                self.allele_variant(value, Self::protein_variant)
            }
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
