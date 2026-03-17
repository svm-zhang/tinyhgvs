use pyo3::prelude::*;
use pyo3::types::{PyModule, PyTuple};
use tinyhgvs::{
    parse_hgvs as parse_hgvs_core, CoordinateSystem, HgvsVariant as CoreHgvsVariant,
    NucleotideEdit, NucleotidePositionAnchor, NucleotideSequence, NucleotideSequenceComponent,
    NucleotideSequenceSegment, ParseHgvsError as CoreParseHgvsError, ParseHgvsErrorKind,
    ProteinEdit, ProteinEffect, ProteinSequence, Range, ReferenceSpec, SequenceId, SequenceKind,
    SequenceSource, VariantDescription,
};

const PY_ERRORS_MODULE: &str = "tinyhgvs.errors";
const PY_MODELS_MODULE: &str = "tinyhgvs.models";

#[pyfunction]
fn parse_hgvs<'py>(py: Python<'py>, input: &str) -> PyResult<Bound<'py, PyAny>> {
    match parse_hgvs_core(input) {
        Ok(variant) => PyModelFactory::import(py)?.variant(&variant),
        Err(error) => Err(PyErrorFactory::import(py)?.parse_error(&error)?),
    }
}

#[pymodule]
fn _tinyhgvs(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_hgvs, m)?)?;
    Ok(())
}

struct PyModelFactory<'py> {
    py: Python<'py>,
    module: Bound<'py, PyModule>,
}

impl<'py> PyModelFactory<'py> {
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

    fn sequence_kind(&self, value: SequenceKind) -> PyResult<Bound<'py, PyAny>> {
        self.class("SequenceKind")?
            .call1((sequence_kind_value(value),))
    }

    fn position_anchor(&self, value: NucleotidePositionAnchor) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotidePositionAnchor")?
            .call1((position_anchor_value(value),))
    }

    fn sequence_id(&self, value: &SequenceId) -> PyResult<Bound<'py, PyAny>> {
        self.class("SequenceId")?.call1((
            &value.raw,
            value.version,
            self.sequence_kind(value.kind)?,
        ))
    }

    fn reference_spec(&self, value: &ReferenceSpec) -> PyResult<Bound<'py, PyAny>> {
        let context = value
            .context
            .as_ref()
            .map(|context| self.sequence_id(context))
            .transpose()?;

        self.class("ReferenceSpec")?
            .call1((self.sequence_id(&value.primary)?, context))
    }

    fn nucleotide_position(
        &self,
        value: &tinyhgvs::NucleotidePosition,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotidePosition")?.call1((
            self.position_anchor(value.anchor)?,
            value.position,
            value.offset,
        ))
    }

    fn protein_position(&self, value: &tinyhgvs::ProteinPosition) -> PyResult<Bound<'py, PyAny>> {
        self.class("ProteinPosition")?
            .call1((&value.residue, value.ordinal))
    }

    fn nucleotide_range(
        &self,
        value: &Range<tinyhgvs::NucleotidePosition>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let end = value
            .end
            .as_ref()
            .map(|end| self.nucleotide_position(end))
            .transpose()?;

        self.class("Range")?
            .call1((self.nucleotide_position(&value.start)?, end))
    }

    fn protein_range(
        &self,
        value: &Range<tinyhgvs::ProteinPosition>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let end = value
            .end
            .as_ref()
            .map(|end| self.protein_position(end))
            .transpose()?;

        self.class("Range")?
            .call1((self.protein_position(&value.start)?, end))
    }

    fn sequence_source(&self, value: &SequenceSource) -> PyResult<Bound<'py, PyAny>> {
        match value {
            SequenceSource::CurrentReference => self.class("CurrentReferenceSource")?.call0(),
            SequenceSource::OtherReference {
                reference,
                coordinate_system,
            } => self.class("OtherReferenceSource")?.call1((
                self.reference_spec(reference)?,
                self.coordinate_system(*coordinate_system)?,
            )),
        }
    }

    fn nucleotide_sequence_segment(
        &self,
        value: &NucleotideSequenceSegment,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotideSequenceSegment")?.call1((
            self.sequence_source(&value.source)?,
            self.nucleotide_range(&value.location)?,
            value.is_inverted,
        ))
    }

    fn nucleotide_sequence_component(
        &self,
        value: &NucleotideSequenceComponent,
    ) -> PyResult<Bound<'py, PyAny>> {
        match value {
            NucleotideSequenceComponent::Literal(literal) => {
                self.class("LiteralSequenceComponent")?.call1((literal,))
            }
            NucleotideSequenceComponent::Repeat { unit, count } => {
                self.class("RepeatSequenceComponent")?.call1((unit, *count))
            }
            NucleotideSequenceComponent::Segment(segment) => self
                .class("SegmentSequenceComponent")?
                .call1((self.nucleotide_sequence_segment(segment)?,)),
        }
    }

    fn nucleotide_sequence(&self, value: &NucleotideSequence) -> PyResult<Bound<'py, PyAny>> {
        let components = value
            .components
            .iter()
            .map(|component| self.nucleotide_sequence_component(component))
            .collect::<PyResult<Vec<_>>>()?;
        let components = PyTuple::new(self.py, components)?;

        self.class("NucleotideSequence")?.call1((components,))
    }

    fn protein_sequence(&self, value: &ProteinSequence) -> PyResult<Bound<'py, PyAny>> {
        let residues = PyTuple::new(self.py, &value.residues)?;
        self.class("ProteinSequence")?.call1((residues,))
    }

    fn nucleotide_edit(&self, value: &NucleotideEdit) -> PyResult<Bound<'py, PyAny>> {
        match value {
            NucleotideEdit::NoChange => self.class("NucleotideNoChangeEdit")?.call0(),
            NucleotideEdit::Substitution {
                reference,
                alternate,
            } => self
                .class("NucleotideSubstitutionEdit")?
                .call1((reference, alternate)),
            NucleotideEdit::Deletion => self.class("NucleotideDeletionEdit")?.call0(),
            NucleotideEdit::Duplication => self.class("NucleotideDuplicationEdit")?.call0(),
            NucleotideEdit::Insertion { sequence } => self
                .class("NucleotideInsertionEdit")?
                .call1((self.nucleotide_sequence(sequence)?,)),
            NucleotideEdit::Inversion => self.class("NucleotideInversionEdit")?.call0(),
            NucleotideEdit::DeletionInsertion { sequence } => self
                .class("NucleotideDeletionInsertionEdit")?
                .call1((self.nucleotide_sequence(sequence)?,)),
        }
    }

    fn protein_edit(&self, value: &ProteinEdit) -> PyResult<Bound<'py, PyAny>> {
        match value {
            ProteinEdit::Unknown => self.class("ProteinUnknownEdit")?.call0(),
            ProteinEdit::NoChange => self.class("ProteinNoChangeEdit")?.call0(),
            ProteinEdit::Substitution { to } => self.class("ProteinSubstitutionEdit")?.call1((to,)),
            ProteinEdit::Deletion => self.class("ProteinDeletionEdit")?.call0(),
            ProteinEdit::Duplication => self.class("ProteinDuplicationEdit")?.call0(),
            ProteinEdit::Insertion { sequence } => self
                .class("ProteinInsertionEdit")?
                .call1((self.protein_sequence(sequence)?,)),
            ProteinEdit::DeletionInsertion { sequence } => self
                .class("ProteinDeletionInsertionEdit")?
                .call1((self.protein_sequence(sequence)?,)),
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
                .call1((self.protein_range(location)?, self.protein_edit(edit)?)),
        }
    }

    fn description(&self, value: &VariantDescription) -> PyResult<Bound<'py, PyAny>> {
        match value {
            VariantDescription::Nucleotide(value) => self.class("NucleotideVariant")?.call1((
                self.nucleotide_range(&value.location)?,
                self.nucleotide_edit(&value.edit)?,
            )),
            VariantDescription::Protein(value) => self
                .class("ProteinVariant")?
                .call1((value.is_predicted, self.protein_effect(&value.effect)?)),
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

fn sequence_kind_value(value: SequenceKind) -> &'static str {
    match value {
        SequenceKind::RefSeqChromosome => "refseq_chromosome",
        SequenceKind::RefSeqContig => "refseq_contig",
        SequenceKind::RefSeqGeneRegion => "refseq_gene_region",
        SequenceKind::RefSeqCodingTranscript => "refseq_coding_transcript",
        SequenceKind::RefSeqNoncodingTranscript => "refseq_noncoding_transcript",
        SequenceKind::RefSeqProtein => "refseq_protein",
        SequenceKind::EnsemblGene => "ensembl_gene",
        SequenceKind::EnsemblTranscript => "ensembl_transcript",
        SequenceKind::EnsemblProtein => "ensembl_protein",
        SequenceKind::LrgGeneRegion => "lrg_gene_region",
        SequenceKind::LrgTranscript => "lrg_transcript",
        SequenceKind::LrgProtein => "lrg_protein",
        SequenceKind::Unknown => "unknown",
    }
}

fn position_anchor_value(value: NucleotidePositionAnchor) -> &'static str {
    match value {
        NucleotidePositionAnchor::Coordinate => "coordinate",
        NucleotidePositionAnchor::CdsStart => "cds_start",
        NucleotidePositionAnchor::CdsEnd => "cds_end",
    }
}

fn error_kind_value(value: ParseHgvsErrorKind) -> &'static str {
    match value {
        ParseHgvsErrorKind::InvalidSyntax => "invalid_syntax",
        ParseHgvsErrorKind::UnsupportedSyntax => "unsupported_syntax",
        ParseHgvsErrorKind::SemanticConstraint => "semantic_constraint",
    }
}
