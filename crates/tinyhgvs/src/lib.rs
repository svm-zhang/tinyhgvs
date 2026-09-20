//! Lightweight HGVS variant parser.
//!
//! `tinyhgvs` parses an HGVS variant into explicit Rust structs and enums
//! describing:
//!
//! - the reference sequence context such as `NM_004006.2` or `NP_003997.1`
//! - the coordinate type such as coding DNA (`c.`), genomic DNA (`g.`), RNA
//!   (`r.`), or protein (`p.`)
//! - the biological description itself, represented as genomic, coding-DNA,
//!   RNA, or protein outcomes and allele forms
//!
//! The main entry points are:
//!
//! - [`parse_hgvs`] to parse a string into [`HgvsVariant`]
//! - [`ParseHgvsError`] to inspect invalid or unsupported input
//!
//! # Reading the Parsed Model
//!
//! The [`HgvsVariant`] separates a HGVS syntax into three top-level parts:
//!
//! - `reference`: the reference source for a variant.
//! - `coordinate_system`: the one-letter HGVS coordinate type.
//! - `description`: the nucleotide or protein variant description, including
//!   location and base edits or effects.
//!
//! # Examples
//!
//! A coding-DNA substitution with an intronic offset:
//!
//! ```rust
//! use tinyhgvs::{
//!     CodingDnaOutcome, CoordinateSystem, NucleotideAnchor, NucleotideEditKind,
//!     VariantDescription, parse_hgvs,
//! };
//!
//! # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
//! let variant = parse_hgvs("NM_004006.2:c.357+1G>A")?;
//! assert_eq!(variant.coordinate_system, CoordinateSystem::CodingDna);
//!
//! let VariantDescription::CodingDna(CodingDnaOutcome::Known(edit)) = variant.description else {
//!     panic!("expected a coding-DNA edit");
//! };
//!
//! assert_eq!(edit.location.start().unwrap().anchor(), Some(NucleotideAnchor::Absolute));
//! assert_eq!(edit.location.start().unwrap().coordinate(), Some(357));
//! assert_eq!(edit.location.start().unwrap().offset(), Some(1));
//! assert!(matches!(
//!     edit.kind,
//!     NucleotideEditKind::Substitution { ref reference, ref alternate }
//!         if reference == "G" && alternate == "A"
//! ));
//! # Ok(())
//! # }
//! ```
//!
//! An RNA special outcome:
//!
//! ```rust
//! use tinyhgvs::{RnaOutcome, VariantDescription, parse_hgvs};
//!
//! # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
//! let variant = parse_hgvs("NM_004006.3:r.spl")?;
//!
//! let VariantDescription::Rna(outcome) = variant.description else {
//!     panic!("expected an RNA outcome");
//! };
//!
//! assert!(matches!(outcome, RnaOutcome::UncertainSplicing));
//! # Ok(())
//! # }
//! ```
//!
//! A nonsense mutation leading to an early termination at protein-level:
//!
//! ```rust
//! use tinyhgvs::{
//!     CoordinateSystem, OutcomeCertainty, ProteinEditKind, ProteinOutcome,
//!     VariantDescription, parse_hgvs,
//! };
//!
//! # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
//! let variant = parse_hgvs("NP_003997.1:p.Trp24Ter")?;
//! assert_eq!(variant.coordinate_system, CoordinateSystem::Protein);
//!
//! let VariantDescription::Protein(ProteinOutcome::Produced { edit, certainty }) =
//!     variant.description
//! else {
//!     panic!("expected a produced protein outcome");
//! };
//!
//! assert_eq!(certainty, OutcomeCertainty::Certain);
//! assert_eq!(edit.location.start().unwrap().residue, "Trp");
//! assert!(matches!(
//!     edit.kind,
//!     ProteinEditKind::Substitution { ref to } if to == "Ter"
//! ));
//! # Ok(())
//! # }
//! ```
//!
//! Unsupported syntax is reported with a stable diagnostic code:
//!
//! ```rust
//! use tinyhgvs::parse_hgvs;
//!
//! let error = parse_hgvs("NC_000023.11:g.pter_qtersup").unwrap_err();
//! assert_eq!(error.code(), "unsupported.telomeric_position");
//! ```

mod diagnostics;
mod error;
mod model;
mod parser;

pub use error::{ParseHgvsError, ParseHgvsErrorKind};
pub use model::{
    Accession, Allele, AlleleForm, AllelePhase, AlleleStateCertainty, AlleleVariant,
    CodingDnaOutcome, CoordinateSystem, CopiedSequenceItem, DerivedAllele, GenomicOutcome,
    HgvsVariant, Interval, LiteralSequenceItem, Location, NucleotideAnchor, NucleotideCoordinate,
    NucleotideEdit, NucleotideEditKind, NucleotideSequenceItem, OutcomeCertainty,
    ProteinCoordinate, ProteinEdit, ProteinEditForm, ProteinEditKind, ProteinExtensionEdit,
    ProteinExtensionTerminal, ProteinFrameshiftStop, ProteinFrameshiftStopKind,
    ProteinInsertionSequence, ProteinOutcome, ProteinSequence, Quantity, ReferenceSpec, RepeatEdit,
    RepeatSequenceUnit, ResidueChange, RnaOutcome, VariantDescription,
};
pub use parser::parse_hgvs;
