//! Core data model for supported HGVS variants.

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
    ProteinCoordinate, ProteinEdit, ProteinEditKind, ProteinExtensionEdit,
    ProteinExtensionTerminal, ProteinFrameshiftStop, ProteinFrameshiftStopKind, ProteinSequence,
};
pub use repeat::{Quantity, RepeatEdit, RepeatSequenceUnit};
