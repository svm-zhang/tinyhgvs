use pyo3::prelude::*;
use pyo3::types::PyTuple;

use tinyhgvs::{
    Accession, AllelePhase, CoordinateSystem, CopiedSequenceItem, Interval, Location,
    NucleotideAnchor, NucleotideCoordinate, NucleotideEdit, NucleotideEditKind,
    NucleotideSequenceItem, ProteinCoordinate, Quantity, ReferenceSpec, RepeatEdit,
    RepeatSequenceUnit,
};

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
            .map(|coordinate_system| self.coordinate_system(*coordinate_system))
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
            self.alleles_tuple(&value.variants_unphased, map_variant)?,
        ))
    }

    fn nucleotide_variant(&self, value: &NucleotideVariant) -> PyResult<Bound<'py, PyAny>> {
        self.class("NucleotideVariant")?.call1((
            self.nucleotide_location(&value.location)?,
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
            ProteinEffect::Known { location, edit } => self
                .class("ProteinEditEffect")?
                .call1((self.protein_location(location)?, self.protein_edit(edit)?)),
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

fn nucleotide_anchor_value(value: NucleotideAnchor) -> &'static str {
    match value {
        NucleotideAnchor::Absolute => "absolute",
        NucleotideAnchor::RelativeCdsStart => "relative_cds_start",
        NucleotideAnchor::RelativeCdsEnd => "relative_cds_end",
    }
}
