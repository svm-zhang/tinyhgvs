//! Core data model for supported HGVS variants.

/// A parsed HGVS variant.
///
/// This is the root model returned by [`crate::parse_hgvs`]. It keeps the
/// reference, coordinate system, and biological description together.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{CoordinateSystem, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_007373.4:c.-1C>T").unwrap();
/// assert_eq!(variant.coordinate_system, CoordinateSystem::CodingDna);
///
/// match variant.description {
///     VariantDescription::Nucleotide(description) => {
///         assert_eq!(description.location.start.coordinate, -1);
///     }
///     VariantDescription::Protein(_) => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HgvsVariant {
    /// Reference sequence metadata such as `NM_004006.2`. Optional for
    /// shorthand protein variants such as `p.Gly2_Met46del`.
    pub reference: Option<ReferenceSpec>,
    /// HGVS coordinate type such as `c`, `g`, `r`, or `p`.
    pub coordinate_system: CoordinateSystem,
    /// Parsed variant description for nucleotide or protein syntax.
    pub description: VariantDescription,
}

/// Reference metadata preceding the `:` in an HGVS expression.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::parse_hgvs;
///
/// let variant = parse_hgvs("NC_000023.11(NM_004006.2):c.3921dup").unwrap();
/// let reference = variant.reference.unwrap();
///
/// assert_eq!(reference.primary.id, "NC_000023.11");
/// assert_eq!(reference.context.unwrap().id, "NM_004006.2");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceSpec {
    /// The main accession being described.
    pub primary: Accession,
    /// Optional contextual accession "NG_012232.1(NM_004006.2):c.93+1G>T"
    pub context: Option<Accession>,
}

/// A parsed accession with optional version.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::Accession;
///
/// let accession = Accession::new("NP_003997.2");
/// assert_eq!(accession.id, "NP_003997.2");
/// assert_eq!(accession.version, Some(2));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Accession {
    /// Accession ID.
    pub id: String,
    /// Optional version for the accession ID.
    pub version: Option<u32>,
}

impl Accession {
    /// Builds an [`Accession`] from an accession string.
    pub fn new(id: impl Into<String>) -> Self {
        let id = id.into();
        let version = id
            .rsplit_once('.')
            .and_then(|(_, suffix)| suffix.parse::<u32>().ok());

        Self { id, version }
    }
}

/// HGVS coordinate system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateSystem {
    Genomic,
    CircularGenomic,
    Mitochondrial,
    CodingDna,
    NonCodingDna,
    Rna,
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
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantDescription {
    Nucleotide(NucleotideVariant),
    Protein(ProteinVariant),
}

/// Parsed nucleotide location and edit.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_004006.2:c.357+1G>A").unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(description) => {
///         assert!(matches!(
///             description.edit,
///             NucleotideEdit::Substitution { ref reference, ref alternate }
///                 if reference == "G" && alternate == "A"
///         ));
///         assert_eq!(description.location.start.coordinate, 357);
///         assert_eq!(description.location.start.offset, 1);
///     }
///     VariantDescription::Protein(_) => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideVariant {
    pub location: Interval<NucleotideCoordinate>,
    pub edit: NucleotideEdit,
}

/// Parsed protein consequence.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.1:p.(Trp24Ter)").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(description) => {
///         assert!(description.is_predicted);
///         assert!(matches!(description.effect, ProteinEffect::Edit { .. }));
///     }
///     VariantDescription::Nucleotide(_) => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinVariant {
    pub is_predicted: bool,
    pub effect: ProteinEffect,
}

/// Supported protein consequence forms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinEffect {
    Unknown,
    NoProteinProduced,
    Edit {
        location: Interval<ProteinCoordinate>,
        edit: ProteinEdit,
    },
}

/// Model describing a stop codon is known (long-form), or omitted (short-form),
/// or unknown (not encountered) due to a frameshift event.
///
/// - "Known" or long-form: `p.Arg97ProfsTer23`
/// - "Omitted" or short-form: `p.Arg97fs`
/// - "Unknown" or "not encountered": `p.Arg97ProfsTer?`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinFrameshiftStopKind {
    Omitted,
    Unknown,
    Known,
}

/// Model describing stop codon information in a protein frameshift edit.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, ProteinFrameshiftStopKind, VariantDescription, parse_hgvs};
///
/// let short = parse_hgvs("NP_0123456.1:p.Arg97fs").unwrap();
/// let known = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23").unwrap();
/// let unknown = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer?").unwrap();
///
/// let extract_stop = |variant: tinyhgvs::HgvsVariant| match variant.description {
///     VariantDescription::Protein(description) => match description.effect {
///         ProteinEffect::Edit { edit: ProteinEdit::Frameshift { stop, .. }, .. } => stop,
///         _ => unreachable!("expected protein frameshift"),
///     },
///     _ => unreachable!("expected protein variant"),
/// };
///
/// let short_stop = extract_stop(short);
/// assert_eq!(short_stop.kind, ProteinFrameshiftStopKind::Omitted);
/// assert_eq!(short_stop.ordinal, None);
///
/// let known_stop = extract_stop(known);
/// assert_eq!(known_stop.kind, ProteinFrameshiftStopKind::Known);
/// assert_eq!(known_stop.ordinal, Some(23));
///
/// let unknown_stop = extract_stop(unknown);
/// assert_eq!(unknown_stop.kind, ProteinFrameshiftStopKind::Unknown);
/// assert_eq!(unknown_stop.ordinal, None);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinFrameshiftStop {
    pub ordinal: Option<usize>,
    pub kind: ProteinFrameshiftStopKind,
}

/// Inclusive interval used for nucleotide and protein locations.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.2:p.Lys23_Val25del").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(description) => {
///         let location = match description.effect {
///             tinyhgvs::ProteinEffect::Edit { ref location, .. } => location,
///             _ => unreachable!("expected protein edit"),
///         };
///         assert_eq!(location.start.residue, "Lys");
///         assert_eq!(location.end.as_ref().unwrap().residue, "Val");
///     }
///     VariantDescription::Nucleotide(_) => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval<T> {
    pub start: T,
    pub end: Option<T>,
}

/// Anchor used by nucleotide coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleotideAnchor {
    /// Coordinate is absolute.
    Absolute,
    /// Coordinate is relative to the CDS start site.
    RelativeCdsStart,
    /// Coordinate is relative to the CDS end site.
    RelativeCdsEnd,
}

/// Nucleotide coordinate with explicit anchor and offset.
///
/// The primary coordinate keeps the sign written in the HGVS string. For
/// example, `c.-1` becomes `coordinate == -1`, while `c.*1` becomes
/// `coordinate == 1`.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{NucleotideAnchor, VariantDescription, parse_hgvs};
///
/// let five_prime = parse_hgvs("NM_007373.4:c.-1C>T").unwrap();
/// let three_prime = parse_hgvs("NM_001272071.2:c.*1C>T").unwrap();
///
/// match five_prime.description {
///     VariantDescription::Nucleotide(description) => {
///         assert_eq!(description.location.start.anchor, NucleotideAnchor::RelativeCdsStart);
///         assert_eq!(description.location.start.coordinate, -1);
///         assert_eq!(description.location.start.offset, 0);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
///
/// match three_prime.description {
///     VariantDescription::Nucleotide(description) => {
///         assert_eq!(description.location.start.anchor, NucleotideAnchor::RelativeCdsEnd);
///         assert_eq!(description.location.start.coordinate, 1);
///         assert_eq!(description.location.start.offset, 0);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideCoordinate {
    /// Anchor type used to describe `coordinate`.
    pub anchor: NucleotideAnchor,
    /// Primary coordinate written in the HGVS string.
    pub coordinate: i32,
    /// Secondary displacement written after the primary coordinate.
    pub offset: i32,
}

impl NucleotideCoordinate {
    /// Returns `true` for splice-adjacent intronic coordinates such as `357+1`.
    pub fn is_intronic(&self) -> bool {
        matches!(self.anchor, NucleotideAnchor::Absolute) && self.offset != 0
    }

    /// Returns `true` for 5' UTR coordinates such as `c.-81`.
    pub fn is_five_prime_utr(&self) -> bool {
        matches!(self.anchor, NucleotideAnchor::RelativeCdsStart)
    }

    /// Returns `true` for 3' UTR coordinates such as `c.*24`.
    pub fn is_three_prime_utr(&self) -> bool {
        matches!(self.anchor, NucleotideAnchor::RelativeCdsEnd)
    }
}

/// Protein coordinate written as amino-acid symbol plus ordinal.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.1:p.Trp24Ter").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(description) => {
///         let location = match description.effect {
///             ProteinEffect::Edit { ref location, .. } => location,
///             _ => unreachable!("expected protein edit"),
///         };
///         assert_eq!(location.start.residue, "Trp");
///         assert_eq!(location.start.ordinal, 24);
///     }
///     _ => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinCoordinate {
    pub residue: String,
    pub ordinal: i32,
}

/// Supported nucleotide edit families.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NucleotideEdit {
    // "="
    NoChange,
    // "G>A"
    Substitution {
        reference: String,
        alternate: String,
    },
    // "del"
    Deletion,
    // "dup"
    Duplication,
    /// Top-level repeated sequence such as `g.123CAG[23]` or
    /// `r.456_465[4]466_489[9]490_499[3]`.
    Repeat {
        blocks: Vec<NucleotideRepeatBlock>,
    },
    Insertion {
        items: Vec<NucleotideSequenceItem>,
    },
    // "inv"
    Inversion,
    DeletionInsertion {
        items: Vec<NucleotideSequenceItem>,
    },
}

/// A single sequence item inside a nucleotide insertion or deletion-insertion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NucleotideSequenceItem {
    Literal(LiteralSequenceItem),
    Repeat(RepeatSequenceItem),
    Copied(CopiedSequenceItem),
}

/// Literal inserted or replacement bases such as `A` or `AGGG`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LiteralSequenceItem {
    pub value: String,
}

/// Repeated inserted or replacement sequence such as `N[12]`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatSequenceItem {
    pub unit: String,
    pub count: usize,
}

/// One repeated block/unit in a nucleotide repeat variant description.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NC_000014.8:g.123CAG[23]").unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(description) => {
///         let NucleotideEdit::Repeat { blocks } = description.edit else {
///             unreachable!("expected nucleotide repeat");
///         };
///         assert_eq!(blocks[0].unit.as_deref(), Some("CAG"));
///         assert_eq!(blocks[0].count, 23);
///     }
///     VariantDescription::Protein(_) => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideRepeatBlock {
    /// Total number of copies reported for the repeated unit.
    pub count: usize,
    /// Explicit repeat unit.
    pub unit: Option<String>,
    /// Repeat unit described as an interval.
    pub location: Option<Interval<NucleotideCoordinate>>,
}

/// Sequence copied from the same or another reference.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, NucleotideSequenceItem, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]").unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(description) => {
///         let NucleotideEdit::Insertion { items } = description.edit else {
///             unreachable!("expected insertion");
///         };
///         let NucleotideSequenceItem::Copied(item) = &items[1] else {
///             unreachable!("expected copied sequence");
///         };
///         assert!(item.is_from_same_reference());
///         assert_eq!(item.source_location.start.coordinate, 450);
///         assert_eq!(item.source_location.end.as_ref().unwrap().coordinate, 470);
///     }
///     VariantDescription::Protein(_) => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CopiedSequenceItem {
    /// `None` when the copied segment comes from the same reference source
    /// as written in the reference metadata.
    pub source_reference: Option<ReferenceSpec>,
    /// `None` when the same coordinate system is used as written in the
    /// reference metadata.
    pub source_coordinate_system: Option<CoordinateSystem>,
    /// Interval on the source reference from which sequence is copied.
    pub source_location: Interval<NucleotideCoordinate>,
    /// Whether the copied sequence is inverted.
    pub is_inverted: bool,
}

impl CopiedSequenceItem {
    /// Returns `true` when the copied sequence comes from the outer reference.
    pub fn is_from_same_reference(&self) -> bool {
        self.source_reference.is_none() && self.source_coordinate_system.is_none()
    }
}

/// Supported protein edit families in the first release.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinEdit {
    Unknown,
    NoChange,
    Substitution {
        to: String,
    },
    Deletion,
    Duplication,
    /// Top-level repeated sequence such as `p.Ala2[10]` or
    /// `p.Arg65_Ser67[12]`.
    Repeat {
        count: usize,
    },
    /// Protein frameshift such as `p.Arg97fs` or `p.Arg97ProfsTer23`.
    Frameshift {
        to_residue: Option<String>,
        stop: ProteinFrameshiftStop,
    },
    Insertion {
        sequence: ProteinSequence,
    },
    DeletionInsertion {
        sequence: ProteinSequence,
    },
}

/// Ordered protein insertion or replacement sequence.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("p.Lys2_Gly3insGlnSerLys").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(description) => {
///         let ProteinEffect::Edit { edit, .. } = description.effect else {
///             unreachable!("expected protein edit");
///         };
///         let ProteinEdit::Insertion { sequence } = edit else {
///             unreachable!("expected protein insertion");
///         };
///         assert_eq!(
///             sequence.residues,
///             vec!["Gln".to_string(), "Ser".to_string(), "Lys".to_string()]
///         );
///     }
///     VariantDescription::Nucleotide(_) => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinSequence {
    pub residues: Vec<String>,
}
