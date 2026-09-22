use pyo3::prelude::*;
use pyo3::types::{PyModule, PyTuple};

use tinyhgvs::{
    Accession, Allele, AlleleForm, AllelePhase, AlleleStateCertainty, AlleleVariant,
    CodingDnaOutcome, CoordinateSystem, CopiedSequenceItem, GenomicOutcome, HgvsVariant, Interval,
    Location, NucleotideAnchor, NucleotideCoordinate, NucleotideEdit, NucleotideEditKind,
    NucleotideSequenceItem, OutcomeCertainty, ProteinCoordinate, ProteinEdit, ProteinEditForm,
    ProteinEditKind, ProteinExtensionTerminal, ProteinFrameshiftStop, ProteinFrameshiftStopKind,
    ProteinInsertionSequence, ProteinOutcome, ProteinSequence, Quantity, ReferenceSpec, RepeatEdit,
    RepeatSequenceUnit, ResidueChange, RnaOutcome, VariantDescription,
};

const PY_MODELS_MODULE: &str = "tinyhgvs.models";

pub(crate) struct PyModelCodec<'py> {
    py: Python<'py>,
    module: Bound<'py, PyModule>,
}

impl<'py> PyModelCodec<'py> {
    pub(crate) fn import(py: Python<'py>) -> PyResult<Self> {
        Ok(Self {
            py,
            module: PyModule::import(py, PY_MODELS_MODULE)?,
        })
    }

    fn class(&self, name: &str) -> PyResult<Bound<'py, PyAny>> {
        self.module.getattr(name)
    }

    fn coordinate_system(&self, value: &CoordinateSystem) -> PyResult<Bound<'py, PyAny>> {
        self.class("CoordinateSystem")?.call1((value.as_str(),))
    }

    fn nucleotide_anchor(&self, value: NucleotideAnchor) -> PyResult<Bound<'py, PyAny>> {
        let value = match value {
            NucleotideAnchor::Absolute => "absolute",
            NucleotideAnchor::RelativeCdsStart => "relative_cds_start",
            NucleotideAnchor::RelativeCdsEnd => "relative_cds_end",
        };

        self.class("NucleotideAnchor")?.call1((value,))
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
        match value {
            NucleotideCoordinate::Known {
                anchor,
                coordinate,
                offset,
            } => self.class("NucleotideCoordinate")?.call1((
                self.nucleotide_anchor(*anchor)?,
                *coordinate,
                *offset,
            )),
            NucleotideCoordinate::Unknown => self.class("NucleotideCoordinate")?.call1((
                self.py.None(),
                self.py.None(),
                self.py.None(),
            )),
        }
    }

    fn protein_coordinate(&self, value: &ProteinCoordinate) -> PyResult<Bound<'py, PyAny>> {
        self.class("ProteinCoordinate")?
            .call1((&value.residue, value.ordinal))
    }

    fn quantity(&self, value: &Quantity) -> PyResult<Bound<'py, PyAny>> {
        match value {
            Quantity::Known { count } => self.class("KnownQuantity")?.call1((*count,)),

            Quantity::Uncertain { lo, hi } => self.class("UncertainQuantity")?.call1((*lo, *hi)),

            Quantity::Unknown => self.class("UnknownQuantity")?.call0(),
        }
    }

    fn possible_range_with<T>(
        &self,
        value: &Interval<T>,
        map_position: fn(&Self, &T) -> PyResult<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let end = value
            .end
            .as_ref()
            .map(|end| map_position(self, end))
            .transpose()?;

        self.class("PossibleRange")?
            .call1((map_position(self, &value.start)?, end))
    }

    fn location_with<T>(
        &self,
        value: &Location<T>,
        map_position: fn(&Self, &T) -> PyResult<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        match value {
            Location::Known(interval) => {
                let end = interval
                    .end
                    .as_ref()
                    .map(|end| map_position(self, end))
                    .transpose()?;
                self.class("KnownLocation")?
                    .call1((map_position(self, &interval.start)?, end))
            }
            Location::Uncertain(interval) => {
                let end = interval
                    .end
                    .as_ref()
                    .map(|range| self.possible_range_with(range, map_position))
                    .transpose()?;
                self.class("UncertainLocation")?.call1((
                    self.possible_range_with(&interval.start, map_position)?,
                    end,
                ))
            }
        }
    }

    fn nucleotide_location(
        &self,
        value: &Location<NucleotideCoordinate>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.location_with(value, Self::nucleotide_coordinate)
    }

    fn known_nucleotide_location_from_interval(
        &self,
        value: &Interval<NucleotideCoordinate>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let start = self.nucleotide_coordinate(&value.start)?;

        let end = value
            .end
            .as_ref()
            .map(|end| self.nucleotide_coordinate(end))
            .transpose()?;

        self.class("KnownLocation")?.call1((start, end))
    }

    fn protein_location(&self, value: &Location<ProteinCoordinate>) -> PyResult<Bound<'py, PyAny>> {
        self.location_with(value, Self::protein_coordinate)
    }

    fn copied_sequence(&self, value: &CopiedSequenceItem) -> PyResult<Bound<'py, PyAny>> {
        let reference = value
            .source_reference
            .as_ref()
            .map(|reference| self.reference_spec(reference))
            .transpose()?;

        let coordinate_system = value
            .source_coordinate_system
            .as_ref()
            .map(|coordinate_system| self.coordinate_system(coordinate_system))
            .transpose()?;

        let location = self.known_nucleotide_location_from_interval(&value.source_location)?;

        self.class("CopiedSequence")?.call1((
            reference,
            coordinate_system,
            location,
            value.is_inverted,
        ))
    }

    fn repeat_sequence_unit(&self, value: &RepeatSequenceUnit) -> PyResult<Bound<'py, PyAny>> {
        match value {
            RepeatSequenceUnit::Known(item) => self.class("KnownRepeatUnit")?.call1((&item.value,)),

            RepeatSequenceUnit::Unknown => self.class("UnknownRepeatUnit")?.call0(),
        }
    }

    fn repeat_edit(&self, value: &RepeatEdit) -> PyResult<Bound<'py, PyAny>> {
        let unit = value
            .unit
            .as_ref()
            .map(|unit| self.repeat_sequence_unit(unit))
            .transpose()?;

        let quantity = self.quantity(&value.quantity)?;

        self.class("Repeat")?.call1((unit, quantity))
    }

    fn nucleotide_sequence_item(
        &self,
        value: &NucleotideSequenceItem,
    ) -> PyResult<Bound<'py, PyAny>> {
        match value {
            NucleotideSequenceItem::Literal(item) => {
                self.class("LiteralSequence")?.call1((&item.value,))
            }

            NucleotideSequenceItem::Repeat(edit) => self.repeat_edit(edit),

            NucleotideSequenceItem::Copied(item) => self.copied_sequence(item),
        }
    }

    fn nucleotide_sequence_items_to_py_tuple(
        &self,
        items: &[NucleotideSequenceItem],
    ) -> PyResult<Bound<'py, PyTuple>> {
        let items = items
            .iter()
            .map(|item| self.nucleotide_sequence_item(item))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(self.py, items)
    }

    fn repeat_blocks_to_py_tuple(&self, blocks: &[RepeatEdit]) -> PyResult<Bound<'py, PyTuple>> {
        let items = blocks
            .iter()
            .map(|item| self.repeat_edit(item))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(self.py, items)
    }

    fn allele_state_certainty(&self, value: &AlleleStateCertainty) -> PyResult<Bound<'py, PyAny>> {
        let value = match value {
            AlleleStateCertainty::Certain => "certain",
            AlleleStateCertainty::Uncertain => "uncertain",
        };

        self.class("AlleleStateCertainty")?.call1((value,))
    }

    fn allele_phase(&self, value: &AllelePhase) -> PyResult<Bound<'py, PyAny>> {
        let value = match value {
            AllelePhase::Trans => "trans",
            AllelePhase::Uncertain => "uncertain",
        };
        self.class("AllelePhase")?.call1((value,))
    }

    fn allele_with<T>(
        &self,
        value: &Allele<T>,
        convert: impl Fn(&Self, &T) -> PyResult<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let variants = value
            .variants
            .iter()
            .map(|variant| convert(self, variant))
            .collect::<PyResult<Vec<_>>>()?;

        let variants = PyTuple::new(self.py, variants)?;
        let state_certainty = self.allele_state_certainty(&value.state_certainty)?;

        self.class("Allele")?.call1((variants, state_certainty))
    }

    fn allele_variant_with<T>(
        &self,
        value: &AlleleVariant<T>,
        convert: impl Copy + Fn(&Self, &T) -> PyResult<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let allele_one = self.allele_with(&value.allele_one, convert)?;

        let allele_two = value
            .allele_two
            .as_ref()
            .map(|allele| self.allele_with(allele, convert))
            .transpose()?;

        let phase = value
            .phase
            .map(|phase| self.allele_phase(&phase))
            .transpose()?;

        let unphased = value
            .variants_unphased
            .iter()
            .map(|variant| convert(self, variant))
            .collect::<PyResult<Vec<_>>>()?;

        let unphased = PyTuple::new(self.py, unphased)?;

        self.class("AlleleVariant")?
            .call1((allele_one, allele_two, phase, unphased))
    }

    fn allele_form_with<T>(
        &self,
        value: &AlleleForm<T>,
        convert: impl Copy + Fn(&Self, &T) -> PyResult<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        match value {
            AlleleForm::Single(variant) => self.allele_variant_with(variant, convert),

            AlleleForm::Derived(derived) => {
                let outcomes = derived
                    .outcomes
                    .iter()
                    .map(|outcome| convert(self, outcome))
                    .collect::<PyResult<Vec<_>>>()?;

                let outcomes = PyTuple::new(self.py, outcomes)?;

                self.class("DerivedAlleleForm")?.call1((outcomes,))
            }

            AlleleForm::Alternative(alternatives) => {
                let alternatives = alternatives
                    .iter()
                    .map(|variant| self.allele_variant_with(variant, convert))
                    .collect::<PyResult<Vec<_>>>()?;

                let alternatives = PyTuple::new(self.py, alternatives)?;

                self.class("AlternativeAlleleForm")?.call1((alternatives,))
            }
        }
    }

    fn genomic_allele_form(
        &self,
        value: &AlleleForm<GenomicOutcome>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.allele_form_with(value, |codec, outcome| match outcome {
            GenomicOutcome::Known(edit) => codec.nucleotide_edit(edit),
        })
    }

    fn coding_dna_allele_form(
        &self,
        value: &AlleleForm<CodingDnaOutcome>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.allele_form_with(value, |codec, outcome| codec.coding_dna_outcome(outcome))
    }

    fn rna_allele_form(&self, value: &AlleleForm<RnaOutcome>) -> PyResult<Bound<'py, PyAny>> {
        self.allele_form_with(value, |codec, outcome| codec.rna_outcome(outcome))
    }

    fn protein_allele_form(
        &self,
        value: &AlleleForm<ProteinOutcome>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.allele_form_with(value, |codec, outcome| codec.protein_outcome(outcome))
    }

    fn residue_change(&self, value: &ResidueChange) -> PyResult<Bound<'py, PyAny>> {
        match value {
            ResidueChange::Known(residue) => Ok(residue.into_pyobject(self.py)?.into_any()),

            ResidueChange::Alternative(residues) => {
                Ok(PyTuple::new(self.py, residues.iter())?.into_any())
            }
        }
    }

    fn protein_sequence_to_py_tuple(
        &self,
        value: &ProteinSequence,
    ) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(self.py, value.residues.iter())
    }

    fn nucleotide_edit(&self, value: &NucleotideEdit) -> PyResult<Bound<'py, PyAny>> {
        let location = self.nucleotide_location(&value.location)?;

        match &value.kind {
            NucleotideEditKind::NoChange => self.class("NucleotideNoChange")?.call1((location,)),

            NucleotideEditKind::Substitution {
                reference,
                alternate,
            } => self
                .class("NucleotideSubstitution")?
                .call1((location, reference, alternate)),

            NucleotideEditKind::Deletion => self.class("NucleotideDeletion")?.call1((location,)),

            NucleotideEditKind::Duplication => {
                self.class("NucleotideDuplication")?.call1((location,))
            }

            NucleotideEditKind::Inversion => self.class("NucleotideInversion")?.call1((location,)),

            NucleotideEditKind::Repeat { blocks } => {
                let sequence = self.repeat_blocks_to_py_tuple(blocks)?;

                self.class("NucleotideRepeat")?.call1((location, sequence))
            }

            NucleotideEditKind::Insertion { items } => {
                let sequence = self.nucleotide_sequence_items_to_py_tuple(items)?;

                self.class("NucleotideInsertion")?
                    .call1((location, sequence))
            }

            NucleotideEditKind::DeletionInsertion { items } => {
                let sequence = self.nucleotide_sequence_items_to_py_tuple(items)?;

                self.class("NucleotideDeletionInsertion")?
                    .call1((location, sequence))
            }
        }
    }

    fn protein_frameshift_stop(
        &self,
        value: &ProteinFrameshiftStop,
    ) -> PyResult<Bound<'py, PyAny>> {
        match value.kind {
            ProteinFrameshiftStopKind::Omitted => {
                self.class("OmittedProteinFrameshiftStop")?.call0()
            }

            ProteinFrameshiftStopKind::Unknown => {
                self.class("UnknownProteinFrameshiftStop")?.call0()
            }

            ProteinFrameshiftStopKind::Known => {
                self.class("KnownProteinFrameshiftStop")?.call1((value
                    .ordinal
                    .expect("known frameshift stop requires ordinal"),))
            }
        }
    }

    fn protein_extension_terminal(
        &self,
        value: ProteinExtensionTerminal,
    ) -> PyResult<Bound<'py, PyAny>> {
        let name = match value {
            ProteinExtensionTerminal::N => "N",
            ProteinExtensionTerminal::C => "C",
        };
        self.class("ProteinExtensionTerminal")?.call1((name,))
    }

    fn protein_edit(&self, value: &ProteinEdit) -> PyResult<Bound<'py, PyAny>> {
        let location = self.protein_location(&value.location)?;

        match &value.kind {
            ProteinEditKind::Substitution { to } => {
                let to = self.residue_change(to)?;

                self.class("ProteinSubstitution")?.call1((location, to))
            }

            ProteinEditKind::Deletion => self.class("ProteinDeletion")?.call1((location,)),

            ProteinEditKind::Duplication => self.class("ProteinDuplication")?.call1((location,)),

            ProteinEditKind::Repeat(repeat) => {
                let repeat = self.repeat_edit(repeat)?;

                self.class("ProteinRepeat")?.call1((location, repeat))
            }

            ProteinEditKind::Extension(extension) => {
                let terminal = self.protein_extension_terminal(extension.to_terminal)?;

                self.class("ProteinExtension")?.call1((
                    location,
                    terminal,
                    extension.to_residue.as_deref(),
                    extension.terminal_ordinal,
                ))
            }

            ProteinEditKind::Frameshift { to_residue, stop } => {
                let to_residue = to_residue
                    .as_ref()
                    .map(|residue| self.residue_change(residue))
                    .transpose()?;

                let stop = self.protein_frameshift_stop(stop)?;

                self.class("ProteinFrameshift")?
                    .call1((location, to_residue, stop))
            }

            ProteinEditKind::Insertion { sequence } => match sequence {
                ProteinInsertionSequence::Known(sequence) => {
                    let sequence = self.protein_sequence_to_py_tuple(sequence)?;

                    self.class("KnownProteinInsertion")?
                        .call1((location, sequence))
                }

                ProteinInsertionSequence::Unknown { count } => self
                    .class("UnknownProteinInsertion")?
                    .call1((location, *count)),

                ProteinInsertionSequence::Terminating { ordinal } => self
                    .class("TerminatingProteinInsertion")?
                    .call1((location, *ordinal)),
            },

            ProteinEditKind::DeletionInsertion { sequence } => {
                let sequence = self.protein_sequence_to_py_tuple(sequence)?;

                self.class("ProteinDeletionInsertion")?
                    .call1((location, sequence))
            }

            ProteinEditKind::NoChange(_) => {
                unreachable!("protein no-change edits must be handled by protein_outcome")
            }
        }
    }

    fn outcome_certainty(&self, value: &OutcomeCertainty) -> PyResult<Bound<'py, PyAny>> {
        let value = match value {
            OutcomeCertainty::Certain => "certain",
            OutcomeCertainty::Predicted => "predicted",
        };

        self.class("OutcomeCertainty")?.call1((value,))
    }

    fn coding_dna_outcome(&self, value: &CodingDnaOutcome) -> PyResult<Bound<'py, PyAny>> {
        match value {
            CodingDnaOutcome::Known(edit) => self.nucleotide_edit(edit),

            CodingDnaOutcome::Unknown => self.class("CodingDnaUnknown")?.call0(),
        }
    }

    fn rna_outcome(&self, value: &RnaOutcome) -> PyResult<Bound<'py, PyAny>> {
        match value {
            RnaOutcome::Produced { edit, certainty } => {
                let edit = self.nucleotide_edit(edit)?;
                let certainty = self.outcome_certainty(certainty)?;

                self.class("RnaProduced")?.call1((edit, certainty))
            }

            RnaOutcome::NoChange(certainty) => {
                let certainty = self.outcome_certainty(certainty)?;

                self.class("RnaNoChange")?.call1((certainty,))
            }

            RnaOutcome::NoneProduced(certainty) => {
                let certainty = self.outcome_certainty(certainty)?;

                self.class("RnaNotProduced")?.call1((certainty,))
            }

            RnaOutcome::UncertainSplicing => self.class("RnaUncertainSplicing")?.call0(),

            RnaOutcome::Unknown => self.class("RnaUnknown")?.call0(),

            RnaOutcome::Indeterminate => self.class("RnaIndeterminate")?.call0(),
        }
    }

    fn protein_outcome(&self, value: &ProteinOutcome) -> PyResult<Bound<'py, PyAny>> {
        match value {
            ProteinOutcome::Unknown => self.class("ProteinUnknown")?.call0(),

            ProteinOutcome::NoneProduced(certainty) => {
                let certainty = self.outcome_certainty(certainty)?;

                self.class("ProteinNotProduced")?.call1((certainty,))
            }

            ProteinOutcome::Produced { edit, certainty } => match edit {
                ProteinEditForm::Single(edit) => match &edit.kind {
                    ProteinEditKind::NoChange(certainty) => {
                        let certainty = self.outcome_certainty(certainty)?;

                        self.class("ProteinNoChange")?.call1((certainty,))
                    }

                    _ => {
                        let edit = self.protein_edit(edit)?;
                        let certainty = self.outcome_certainty(certainty)?;

                        self.class("ProteinProduced")?.call1((edit, certainty))
                    }
                },

                ProteinEditForm::Alternative(edits) => {
                    let edits = edits
                        .iter()
                        .map(|edit| self.protein_edit(edit))
                        .collect::<PyResult<Vec<_>>>()?;

                    let edits = PyTuple::new(self.py, edits)?;
                    let certainty = self.outcome_certainty(certainty)?;

                    self.class("ProteinProducedAlternatives")?
                        .call1((edits, certainty))
                }
            },
        }
    }

    fn variant_description(&self, value: &VariantDescription) -> PyResult<Bound<'py, PyAny>> {
        match value {
            VariantDescription::Genomic(outcome) => match outcome {
                GenomicOutcome::Known(edit) => self.nucleotide_edit(edit),
            },

            VariantDescription::CodingDna(outcome) => {
                // Placeholder:
                // use the CURRENT existing cDNA conversion.
                self.coding_dna_outcome(outcome)
            }

            VariantDescription::Rna(outcome) => self.rna_outcome(outcome),

            VariantDescription::Protein(outcome) => self.protein_outcome(outcome),

            VariantDescription::GenomicAllele(form) => self.genomic_allele_form(form),

            VariantDescription::CodingDnaAllele(form) => self.coding_dna_allele_form(form),

            VariantDescription::RnaAllele(form) => self.rna_allele_form(form),

            VariantDescription::ProteinAllele(form) => self.protein_allele_form(form),
        }
    }

    pub(crate) fn hgvs_variant(&self, value: &HgvsVariant) -> PyResult<Bound<'py, PyAny>> {
        let reference = value
            .reference
            .as_ref()
            .map(|reference| self.reference_spec(reference))
            .transpose()?;

        let coordinate_system = self.coordinate_system(&value.coordinate_system)?;

        let description = self.variant_description(&value.description)?;

        self.class("HgvsVariant")?
            .call1((reference, coordinate_system, description))
    }
}
