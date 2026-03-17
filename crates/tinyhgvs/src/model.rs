//! Core data model for supported HGVS variants.
//!
//! The model separates an HGVS expression into biologically meaningful layers:
//!
//! - reference metadata such as [`ReferenceSpec`] and [`SequenceId`]
//! - the coordinate system such as genomic DNA, coding DNA, RNA, or protein
//! - positions and ranges that describe where the change occurs
//! - edit descriptions that describe what sequence change or protein consequence
//!   is being asserted
//!
//! The current design stays intentionally compact. It represents the subset
//! already supported by the parser, while leaving room for the unsupported
//! families cataloged in Phase 4.
//!
//! # Example
//!
//! The coding-DNA variant `NM_004006.2:c.357+1G>A` parses into:
//!
//! - a reference accession `NM_004006.2`
//! - a coding-DNA coordinate system
//! - a nucleotide location anchored at coordinate `357` with intronic offset
//!   `+1`
//! - a substitution from `G` to `A`
//!
//! ```rust
//! use tinyhgvs::{NucleotideEdit, NucleotidePositionAnchor, VariantDescription, parse_hgvs};
//!
//! let variant = parse_hgvs("NM_004006.2:c.357+1G>A").unwrap();
//!
//! match variant.description {
//!     VariantDescription::Nucleotide(nucleotide) => {
//!         assert_eq!(nucleotide.location.start.anchor, NucleotidePositionAnchor::Coordinate);
//!         assert_eq!(nucleotide.location.start.position, Some(357));
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

/// A parsed HGVS variant.
///
/// This is the top-level Rust representation returned by [`crate::parse_hgvs`].
/// It combines:
///
/// - the optional reference sequence prefix
/// - the HGVS coordinate system
/// - the biological description of the variant itself
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{CoordinateSystem, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.1:p.Trp24Ter").unwrap();
/// assert_eq!(variant.coordinate_system, CoordinateSystem::Protein);
/// assert!(variant.reference.is_some());
/// assert!(matches!(variant.description, VariantDescription::Protein(_)));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HgvsVariant {
    /// Optional reference sequence prefix such as `NM_004006.2`.
    pub reference: Option<ReferenceSpec>,
    /// HGVS coordinate system marker such as `c`, `g`, `r`, or `p`.
    pub coordinate_system: CoordinateSystem,
    /// Parsed variant description for nucleotide or protein syntax.
    pub description: VariantDescription,
}

/// Reference metadata preceding the `:` in an HGVS expression.
///
/// This corresponds to the accession part before the coordinate system. It may
/// be a single accession, such as `NM_004006.2`, or a primary accession plus a
/// contextual accession in parentheses, such as
/// `NC_000023.11(NM_004006.2)`.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::parse_hgvs;
///
/// let variant = parse_hgvs("NC_000023.11(NM_004006.2):g.32332592G>T").unwrap();
/// let reference = variant.reference.unwrap();
///
/// assert_eq!(reference.primary.raw, "NC_000023.11");
/// assert_eq!(reference.context.unwrap().raw, "NM_004006.2");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceSpec {
    /// The main accession being described.
    pub primary: SequenceId,
    /// Optional contextual accession such as a transcript inside genomic notation.
    pub context: Option<SequenceId>,
}

/// A parsed accession with inferred sequence kind and optional version.
///
/// A `SequenceId` keeps the accession exactly as written in the HGVS string,
/// while also extracting two useful pieces of metadata:
///
/// - the terminal numeric version, when present
/// - the broad accession family, such as RefSeq coding transcript or Ensembl
///   protein
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{SequenceId, SequenceKind};
///
/// let transcript = SequenceId::new("NM_004006.2");
/// assert_eq!(transcript.raw, "NM_004006.2");
/// assert_eq!(transcript.version, Some(2));
/// assert_eq!(transcript.kind, SequenceKind::RefSeqCodingTranscript);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceId {
    /// Original accession string as written by the user.
    pub raw: String,
    /// Parsed numeric suffix after the final `.` when present.
    pub version: Option<u32>,
    /// Inferred accession family used for downstream dispatch and display.
    pub kind: SequenceKind,
}

impl SequenceId {
    /// Builds a [`SequenceId`] from a raw accession string.
    ///
    /// The accession family is inferred heuristically from the prefix.
    pub fn new(raw: impl Into<String>) -> Self {
        let raw = raw.into();
        let version = raw
            .rsplit_once('.')
            .and_then(|(_, suffix)| suffix.parse::<u32>().ok());

        Self {
            kind: SequenceKind::from_accession(&raw),
            version,
            raw,
        }
    }
}

/// Known accession families inferred from HGVS references.
///
/// These families are intentionally coarse. They help explain what kind of
/// sequence an accession refers to without committing to deeper database
/// semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceKind {
    /// RefSeq chromosome such as `NC_000023.11`.
    RefSeqChromosome,
    /// RefSeq contig such as `NT_` or `NW_`.
    RefSeqContig,
    /// RefSeq gene region such as `NG_`.
    RefSeqGeneRegion,
    /// RefSeq coding transcript such as `NM_`.
    RefSeqCodingTranscript,
    /// RefSeq non-coding transcript such as `NR_`.
    RefSeqNoncodingTranscript,
    /// RefSeq protein accession such as `NP_`.
    RefSeqProtein,
    /// Ensembl gene accession such as `ENSG...`.
    EnsemblGene,
    /// Ensembl transcript accession such as `ENST...`.
    EnsemblTranscript,
    /// Ensembl protein accession such as `ENSP...`.
    EnsemblProtein,
    /// LRG genomic accession.
    LrgGeneRegion,
    /// LRG transcript accession.
    LrgTranscript,
    /// LRG protein accession.
    LrgProtein,
    /// Accession family not yet recognized.
    Unknown,
}

impl SequenceKind {
    fn from_accession(accession: &str) -> Self {
        if accession.starts_with("NC_") {
            Self::RefSeqChromosome
        } else if accession.starts_with("NT_") || accession.starts_with("NW_") {
            Self::RefSeqContig
        } else if accession.starts_with("NG_") {
            Self::RefSeqGeneRegion
        } else if accession.starts_with("NM_") {
            Self::RefSeqCodingTranscript
        } else if accession.starts_with("NR_") {
            Self::RefSeqNoncodingTranscript
        } else if accession.starts_with("NP_") {
            Self::RefSeqProtein
        } else if accession.starts_with("ENSG") {
            Self::EnsemblGene
        } else if accession.starts_with("ENST") {
            Self::EnsemblTranscript
        } else if accession.starts_with("ENSP") {
            Self::EnsemblProtein
        } else if accession.starts_with("LRG_") && accession.contains('t') {
            Self::LrgTranscript
        } else if accession.starts_with("LRG_") && accession.contains('p') {
            Self::LrgProtein
        } else if accession.starts_with("LRG_") {
            Self::LrgGeneRegion
        } else {
            Self::Unknown
        }
    }
}

/// Supported HGVS coordinate systems.
///
/// This enum describes the sequence layer on which the variant is written.
/// The same biological event may be represented differently at the genomic,
/// transcript, RNA, or protein level.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{CoordinateSystem, parse_hgvs};
///
/// let variant = parse_hgvs("NM_004006.2:c.76A>T").unwrap();
/// assert_eq!(variant.coordinate_system, CoordinateSystem::CodingDna);
/// assert_eq!(variant.coordinate_system.as_str(), "c");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateSystem {
    /// Genomic DNA (`g.`).
    Genomic,
    /// Circular genomic DNA (`o.`).
    CircularGenomic,
    /// Mitochondrial DNA (`m.`).
    Mitochondrial,
    /// Coding DNA (`c.`).
    CodingDna,
    /// Non-coding DNA (`n.`).
    NonCodingDna,
    /// RNA (`r.`).
    Rna,
    /// Protein (`p.`).
    Protein,
}

impl CoordinateSystem {
    /// Returns `true` when the coordinate system is protein-based.
    pub fn is_protein(self) -> bool {
        matches!(self, Self::Protein)
    }

    /// Returns the one-letter HGVS coordinate marker.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Genomic => "g",
            Self::CircularGenomic => "o",
            Self::Mitochondrial => "m",
            Self::CodingDna => "c",
            Self::NonCodingDna => "n",
            Self::Rna => "r",
            Self::Protein => "p",
        }
    }
}

/// Top-level variant description for nucleotide or protein syntax.
///
/// This is the branch point between nucleotide-like descriptions such as
/// `c.76A>T` and protein consequences such as `p.Trp24Ter`.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{VariantDescription, parse_hgvs};
///
/// let coding = parse_hgvs("NM_004006.2:c.76A>T").unwrap();
/// assert!(matches!(coding.description, VariantDescription::Nucleotide(_)));
///
/// let protein = parse_hgvs("NP_003997.1:p.Trp24Ter").unwrap();
/// assert!(matches!(protein.description, VariantDescription::Protein(_)));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantDescription {
    /// A nucleotide-based description such as `c.76A>T`.
    Nucleotide(NucleotideVariant),
    /// A protein-based description such as `p.(Trp24Ter)`.
    Protein(ProteinVariant),
}

/// Parsed nucleotide location and edit.
///
/// This represents a nucleotide-level HGVS description after the parser has
/// separated where the change occurs from what the sequence change is.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_004006.2:c.76_78del").unwrap();
/// let description = variant.description;
///
/// match description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         assert_eq!(nucleotide.location.start.position, Some(76));
///         assert_eq!(nucleotide.location.end.unwrap().position, Some(78));
///         assert_eq!(nucleotide.edit, NucleotideEdit::Deletion);
///     }
///     VariantDescription::Protein(_) => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideVariant {
    /// Edited nucleotide interval.
    pub location: Range<NucleotidePosition>,
    /// Edit applied at `location`.
    pub edit: NucleotideEdit,
}

/// Parsed protein consequence.
///
/// Protein HGVS descriptions can be observed directly or predicted from a
/// nucleotide change. The `is_predicted` flag records whether the original HGVS
/// expression wrapped the consequence in parentheses.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("p.(Trp24Ter)").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(protein) => {
///         assert!(protein.is_predicted);
///         assert!(matches!(protein.effect, ProteinEffect::Edit { .. }));
///     }
///     VariantDescription::Nucleotide(_) => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinVariant {
    /// Whether the protein consequence is predicted and therefore wrapped in `()`.
    pub is_predicted: bool,
    /// Protein effect description.
    pub effect: ProteinEffect,
}

/// Supported protein consequence forms.
///
/// The current first-release protein model supports:
///
/// - unknown protein consequence (`p.?`)
/// - no protein produced (`p.0`)
/// - a concrete amino-acid edit at a known location
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.1:p.Trp24Ter").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(protein) => {
///         assert!(matches!(protein.effect, ProteinEffect::Edit { .. }));
///     }
///     VariantDescription::Nucleotide(_) => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinEffect {
    /// `p.?`
    Unknown,
    /// `p.0`
    NoProteinProduced,
    /// Concrete protein location plus edit such as `p.Trp24Ter`.
    Edit {
        /// Amino-acid interval affected by the edit.
        location: Range<ProteinPosition>,
        /// Protein edit at the given location.
        edit: ProteinEdit,
    },
}

/// Generic start/end range used for nucleotide and protein locations.
///
/// HGVS intervals are inclusive. A single-position description has `end =
/// None`, while an interval such as `76_78` or `Trp24_Cys26` has both a start
/// and an end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Range<T> {
    /// First position in the range.
    pub start: T,
    /// Optional inclusive end position for interval syntax.
    pub end: Option<T>,
}

/// Anchor used by nucleotide positions in coding and transcript-relative syntax.
///
/// These anchors explain how to interpret a numeric value in coding-DNA HGVS:
///
/// - coordinate-relative positions such as `357` or `357+1`
/// - positions upstream of the coding sequence start such as `-18`
/// - positions downstream of the coding sequence end such as `*24`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleotidePositionAnchor {
    /// A direct coordinate such as `357` or `357+1`.
    Coordinate,
    /// Position relative to the CDS start such as `-18`.
    CdsStart,
    /// Position relative to the CDS end such as `*18`.
    CdsEnd,
}

/// Nucleotide position with anchor and optional intronic offset.
///
/// This model covers three common biological situations:
///
/// - a direct nucleotide coordinate such as `76`
/// - a splice-adjacent intronic position such as `357+1` or `357-2`
/// - a coding-DNA UTR-relative position such as `-81` or `*24`
///
/// The parser normalizes upstream coding-DNA positions like `c.-81` to
/// `anchor = CdsStart`, `position = Some(0)`, and `offset = -81`.
/// Downstream coding-DNA positions like `c.*24` become
/// `anchor = CdsEnd`, `position = None`, and `offset = 24`.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{NucleotidePositionAnchor, VariantDescription, parse_hgvs};
///
/// let intronic_variant = parse_hgvs("NM_004006.2:c.357+1G>A").unwrap();
/// let utr_variant = parse_hgvs("NM_004006.2:c.*24A>G").unwrap();
///
/// let intronic_position = match intronic_variant.description {
///     VariantDescription::Nucleotide(nucleotide) => nucleotide.location.start,
///     VariantDescription::Protein(_) => unreachable!("expected nucleotide variant"),
/// };
/// assert_eq!(intronic_position.anchor, NucleotidePositionAnchor::Coordinate);
/// assert_eq!(intronic_position.position, Some(357));
/// assert_eq!(intronic_position.offset, 1);
/// assert!(intronic_position.is_intronic());
///
/// let utr_position = match utr_variant.description {
///     VariantDescription::Nucleotide(nucleotide) => nucleotide.location.start,
///     VariantDescription::Protein(_) => unreachable!("expected nucleotide variant"),
/// };
/// assert_eq!(utr_position.anchor, NucleotidePositionAnchor::CdsEnd);
/// assert_eq!(utr_position.position, None);
/// assert_eq!(utr_position.offset, 24);
/// assert!(utr_position.is_three_prime_utr());
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotidePosition {
    /// Anchor used to interpret `position`.
    pub anchor: NucleotidePositionAnchor,
    /// Base coordinate relative to `anchor` when it is explicitly known.
    ///
    /// For upstream coding-DNA positions such as `c.-81`, the anchor is the
    /// CDS start site and `position` is normalized to `Some(0)` with
    /// `offset = -81`.
    ///
    /// For downstream coding-DNA positions such as `c.*24`, the true terminal
    /// CDS coordinate is not recoverable from the string alone, so `position`
    /// is `None` and `offset = 24`.
    pub position: Option<i32>,
    /// Optional signed offset relative to `position` or the anchor.
    pub offset: i32,
}

impl NucleotidePosition {
    /// Returns `true` for coordinate-anchored positions with a non-zero offset.
    ///
    /// These are the splice-adjacent intronic positions such as `357+1` and
    /// `357-2`.
    pub fn is_intronic(&self) -> bool {
        matches!(self.anchor, NucleotidePositionAnchor::Coordinate) && self.offset != 0
    }

    /// Returns `true` for 5' UTR positions anchored to the CDS start site.
    pub fn is_five_prime_utr(&self) -> bool {
        matches!(self.anchor, NucleotidePositionAnchor::CdsStart)
    }

    /// Returns `true` for 3' UTR positions anchored to the CDS end site.
    pub fn is_three_prime_utr(&self) -> bool {
        matches!(self.anchor, NucleotidePositionAnchor::CdsEnd)
    }
}

/// Protein position written as amino-acid symbol plus ordinal.
///
/// Protein HGVS positions pair an amino-acid symbol with its ordinal, such as
/// `Trp24` or `Gly12`.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.1:p.Trp24Ter").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(protein) => {
///         match protein.effect {
///             ProteinEffect::Edit { location, .. } => {
///                 assert_eq!(location.start.residue, "Trp");
///                 assert_eq!(location.start.ordinal, 24);
///             }
///             _ => unreachable!("expected protein edit"),
///         }
///     }
///     VariantDescription::Nucleotide(_) => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinPosition {
    /// Amino-acid symbol as written in the HGVS input.
    pub residue: String,
    /// Amino-acid ordinal.
    pub ordinal: i32,
}

/// Supported nucleotide edit families.
///
/// The nucleotide edit describes the sequence change at a nucleotide location.
/// It is independent from whether the surrounding coordinate system is genomic,
/// transcript, coding-DNA, or RNA.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, VariantDescription, parse_hgvs};
///
/// let substitution = parse_hgvs("NM_004006.2:c.76A>T").unwrap();
/// let insertion = parse_hgvs("NM_004006.2:c.76_77insA").unwrap();
///
/// match substitution.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         assert!(matches!(nucleotide.edit, NucleotideEdit::Substitution { .. }));
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
///
/// match insertion.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         assert!(matches!(nucleotide.edit, NucleotideEdit::Insertion { .. }));
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NucleotideEdit {
    /// No change (`=`).
    NoChange,
    /// Substitution such as `A>T`.
    Substitution {
        /// Reference residue or bases.
        reference: String,
        /// Alternate residue or bases.
        alternate: String,
    },
    /// Deletion (`del`).
    Deletion,
    /// Duplication (`dup`).
    Duplication,
    /// Insertion (`ins...`).
    Insertion {
        /// Inserted sequence description.
        sequence: NucleotideSequence,
    },
    /// Inversion (`inv`).
    Inversion,
    /// Replacement after deletion (`delins...`).
    DeletionInsertion {
        /// Replacement sequence description.
        sequence: NucleotideSequence,
    },
}

/// Ordered inserted or replacement nucleotide sequence.
///
/// A nucleotide insertion or `delins` sequence may contain one or more
/// components written in order. Those components can be:
///
/// - literal bases
/// - repeat notation such as `N[12]`
/// - copied segments from the current or another reference
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_004006.2:c.76_77ins[AT;N[12]]").unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         match nucleotide.edit {
///             NucleotideEdit::Insertion { sequence } => {
///                 assert_eq!(sequence.components.len(), 2);
///             }
///             _ => unreachable!("expected insertion"),
///         }
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideSequence {
    /// Sequence components written inside the inserted or replacement sequence.
    pub components: Vec<NucleotideSequenceComponent>,
}

/// Supported sequence components for nucleotide insertions and delins.
///
/// This enum lets the model distinguish between literal inserted sequence,
/// repeated sequence units, and copied segments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NucleotideSequenceComponent {
    /// Literal bases such as `A` or `ATC`.
    Literal(String),
    /// A repeat unit such as `T[12]`.
    Repeat { unit: String, count: usize },
    /// A local or remote segment reference.
    Segment(NucleotideSequenceSegment),
}

/// Segment copied from the current or another reference.
///
/// This describes inserted sequence segments such as:
///
/// - `850_900` copied from the current reference
/// - `850_900inv` copied from the current reference in inverted orientation
/// - `NC_000022.10:g.35788169_35788352` copied from another accession
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, NucleotideSequenceComponent, SequenceSource, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs(
///     "NM_004006.2:c.76_77ins[NC_000022.10:g.35788169_35788352;850_900inv]"
/// ).unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         match nucleotide.edit {
///             NucleotideEdit::Insertion { sequence } => {
///                 assert!(matches!(
///                     &sequence.components[0],
///                     NucleotideSequenceComponent::Segment(segment)
///                         if matches!(segment.source, SequenceSource::OtherReference { .. })
///                 ));
///                 assert!(matches!(
///                     &sequence.components[1],
///                     NucleotideSequenceComponent::Segment(segment)
///                         if matches!(segment.source, SequenceSource::CurrentReference)
///                             && segment.is_inverted
///                 ));
///             }
///             _ => unreachable!("expected insertion"),
///         }
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideSequenceSegment {
    /// Source reference of the copied segment.
    pub source: SequenceSource,
    /// Location of the copied segment.
    pub location: Range<NucleotidePosition>,
    /// Whether the copied segment is inverted.
    pub is_inverted: bool,
}

/// Source of a sequence segment embedded inside a nucleotide sequence.
///
/// Segment source is either:
///
/// - the current HGVS reference
/// - another accession plus its coordinate system
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SequenceSource {
    /// Segment copied from the current HGVS reference.
    CurrentReference,
    /// Segment copied from another accession and coordinate system.
    OtherReference {
        /// External reference specification.
        reference: ReferenceSpec,
        /// Coordinate system used by the external reference.
        coordinate_system: CoordinateSystem,
    },
}

/// Supported protein edit families in the first release.
///
/// These edits describe the amino-acid consequence at the protein layer.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, VariantDescription, parse_hgvs};
///
/// let substitution = parse_hgvs("NP_003997.1:p.Trp24Ter").unwrap();
/// let insertion = parse_hgvs("NP_003997.1:p.Lys2_Gly3insSerThr").unwrap();
///
/// match substitution.description {
///     VariantDescription::Protein(protein) => match protein.effect {
///         ProteinEffect::Edit { edit, .. } => {
///             assert!(matches!(edit, ProteinEdit::Substitution { ref to } if to == "Ter"));
///         }
///         _ => unreachable!("expected protein edit"),
///     },
///     _ => unreachable!("expected protein variant"),
/// }
///
/// match insertion.description {
///     VariantDescription::Protein(protein) => match protein.effect {
///         ProteinEffect::Edit { edit, .. } => {
///             assert!(matches!(edit, ProteinEdit::Insertion { .. }));
///         }
///         _ => unreachable!("expected protein edit"),
///     },
///     _ => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinEdit {
    /// Unknown edit (`?`) at a known protein location.
    Unknown,
    /// No change (`=`).
    NoChange,
    /// Substitution to a different residue symbol.
    Substitution { to: String },
    /// Deletion (`del`).
    Deletion,
    /// Duplication (`dup`).
    Duplication,
    /// Insertion of one or more residues.
    Insertion { sequence: ProteinSequence },
    /// Replacement after deletion.
    DeletionInsertion { sequence: ProteinSequence },
}

/// Ordered protein insertion or replacement sequence.
///
/// Protein insertions and `delins` consequences store the inserted residues in
/// the order written in the HGVS string.
///
/// # Example
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.1:p.Lys2_Gly3insSerThr").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(protein) => match protein.effect {
///         ProteinEffect::Edit { edit, .. } => match edit {
///             ProteinEdit::Insertion { sequence } => {
///                 assert_eq!(sequence.residues, vec!["Ser".to_string(), "Thr".to_string()]);
///             }
///             _ => unreachable!("expected insertion"),
///         },
///         _ => unreachable!("expected protein edit"),
///     },
///     _ => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinSequence {
    /// Amino-acid symbols in the order written.
    pub residues: Vec<String>,
}
