//! Public data model.
//!
//! The submodules under `model/` are private organization details. Users should
//! import model types from the crate root, for example `tinyhgvs::HgvsVariant`
//! and `tinyhgvs::VariantDescription`.

mod allele;
mod core;
mod location;
mod nucleotide;
mod protein;
mod repeat;

pub use allele::{
    Allele, AlleleForm, AllelePhase, AlleleStateCertainty, AlleleVariant, DerivedAllele,
};
pub use core::{
    Accession, CodingDnaOutcome, CoordinateSystem, GenomicOutcome, HgvsVariant, OutcomeCertainty,
    ProteinOutcome, ReferenceSpec, RnaOutcome, VariantDescription,
};
pub use location::{Interval, Location};
pub use nucleotide::{
    CopiedSequenceItem, LiteralSequenceItem, NucleotideAnchor, NucleotideCoordinate,
    NucleotideEdit, NucleotideEditKind, NucleotideSequenceItem,
};
pub use protein::{
    ProteinCoordinate, ProteinEdit, ProteinEditForm, ProteinEditKind, ProteinExtensionEdit,
    ProteinExtensionTerminal, ProteinFrameshiftStop, ProteinFrameshiftStopKind,
    ProteinInsertionSequence, ProteinSequence, ResidueChange,
};
pub use repeat::{Quantity, RepeatEdit, RepeatSequenceUnit};
