//! Lightweight HGVS variant parser.
//!
//! `tinyhgvs` parses a HGVS variant into explicit Rust structs and enums
//! that describe:
//!
//! - the reference sequence context such as `NM_004006.2` or `NP_003997.1`
//! - the coordinate type such as coding DNA (`c.`), genomic DNA (`g.`), RNA
//!   (`r.`), or protein (`p.`)
//! - the biological description itself, represented as either a nucleotide
//!   variant or a protein consequence
//!
//! The crate is intentionally small. It aims to represent common, high-value
//! HGVS syntax clearly, while returning structured errors for syntax families
//! tracked in the unsupported inventory.
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
//!     location and base edits or effects.
//!
//! # Examples
//!
//! A substitution crossing exon/intron border (intronic):
//!
//! ```rust
//! use tinyhgvs::{NucleotideAnchor, NucleotideEdit, VariantDescription, parse_hgvs};
//!
//! let variant = parse_hgvs("NM_004006.2:c.357+1G>A").unwrap();
//! let description = variant.description;
//!
//! match description {
//!     VariantDescription::Nucleotide(nucleotide) => {
//!         assert_eq!(nucleotide.location.start.anchor, NucleotideAnchor::Absolute);
//!         assert_eq!(nucleotide.location.start.coordinate, 357);
//!         assert_eq!(nucleotide.location.start.offset, 1);
//!         assert!(matches!(
//!             nucleotide.edit,
//!             NucleotideEdit::Substitution { ref reference, ref alternate }
//!                 if reference == "G" && alternate == "A"
//!         ));
//!     }
//!     VariantDescription::Protein(_) => unreachable!("expected nucleotide variant"),
//! }
//! ```
//!
//! A nonsense mutation leading to an early termination at protein-level:
//!
//! ```rust
//! use tinyhgvs::{CoordinateSystem, ProteinEffect, VariantDescription, parse_hgvs};
//!
//! let variant = parse_hgvs("NP_003997.1:p.Trp24Ter").unwrap();
//! assert_eq!(variant.coordinate_system, CoordinateSystem::Protein);
//!
//! match variant.description {
//!     VariantDescription::Protein(protein) => {
//!         assert!(!protein.is_predicted);
//!         assert!(matches!(protein.effect, ProteinEffect::Edit { .. }));
//!     }
//!     VariantDescription::Nucleotide(_) => unreachable!("expected protein variant"),
//! }
//! ```
//!
//! Unsupported syntax is reported with a stable diagnostic code:
//!
//! ```rust
//! use tinyhgvs::parse_hgvs;
//!
//! let error = parse_hgvs("NC_000001.11:g.[123G>A;345del]").unwrap_err();
//! assert_eq!(error.code(), "unsupported.allele");
//! ```

mod diagnostics;
mod error;
mod model;
mod parser;

pub use error::{ParseHgvsError, ParseHgvsErrorKind};
pub use model::{
    Accession, CoordinateSystem, CopiedSequenceItem, HgvsVariant, Interval, LiteralSequenceItem,
    NucleotideAnchor, NucleotideCoordinate, NucleotideEdit, NucleotideSequenceItem,
    NucleotideVariant, ProteinCoordinate, ProteinEdit, ProteinEffect, ProteinSequence,
    ProteinVariant, ReferenceSpec, RepeatSequenceItem, VariantDescription,
};
pub use parser::parse_hgvs;
