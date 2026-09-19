//! Nucleotide coordinates, edits, and inserted sequence items.

use super::{CoordinateSystem, Interval, Location, ReferenceSpec, RepeatEdit};

/// A nucleotide edit applied at a location.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideEdit {
    pub location: Location<NucleotideCoordinate>,
    pub kind: NucleotideEditKind,
}

/// Supported nucleotide edit families.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NucleotideEditKind {
    // c.2376=
    NoChange,
    // c.357+1G>A
    Substitution {
        reference: String,
        alternate: String,
    },
    // c.4072_5145del
    Deletion,
    // g.1234_2345dup
    Duplication,
    /// Top-level repeated sequence such as `g.123CAG[23]`
    Repeat {
        blocks: Vec<RepeatEdit>,
    },
    // c.419_420ins[T;450_470;AGGG]
    Insertion {
        items: Vec<NucleotideSequenceItem>,
    },
    // g.32361330_32361333inv
    Inversion,
    // c.812_829delinsN[12]
    DeletionInsertion {
        items: Vec<NucleotideSequenceItem>,
    },
}

/// A single sequence item inside a nucleotide insertion or deletion-insertion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NucleotideSequenceItem {
    // AGGG
    Literal(LiteralSequenceItem),
    // N[12], CAG[23]
    Repeat(RepeatEdit),
    // 450_470, NC_000022.10:g.35788169_35788352
    Copied(CopiedSequenceItem),
}

/// Literal inserted or replacement bases such as `A` or `AGGG`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LiteralSequenceItem {
    pub value: String,
}

/// Sequence copied from the same or another reference.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{
///     CodingDnaOutcome, NucleotideEditKind, NucleotideSequenceItem, VariantDescription,
///     parse_hgvs,
/// };
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")?;
///
/// let VariantDescription::CodingDna(CodingDnaOutcome::Known(edit)) = variant.description else {
///     panic!("expected a coding-DNA edit");
/// };
///
/// let NucleotideEditKind::Insertion { items } = edit.kind else {
///     panic!("expected insertion");
/// };
/// let NucleotideSequenceItem::Copied(item) = &items[1] else {
///     panic!("expected copied sequence");
/// };
///
/// assert!(item.is_from_same_reference());
/// assert_eq!(item.source_location.start.coordinate(), Some(450));
/// assert_eq!(item.source_location.end.as_ref().unwrap().coordinate(), Some(470));
/// # Ok(())
/// # }
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

/// Anchor used by nucleotide coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleotideAnchor {
    /// Coordinate is absolute, such as `g.123` or `c.357+1`.
    Absolute,
    /// Coordinate is relative to the CDS start site, such as `c.-1`.
    RelativeCdsStart,
    /// Coordinate is relative to the CDS end site, such as `c.*1`.
    RelativeCdsEnd,
}

/// Nucleotide coordinate written as a known position or `?`.
///
/// Known coordinates keep the sign written in the HGVS string. For example,
/// `c.-1` becomes `coordinate() == Some(-1)`, while `c.*1` becomes
/// `coordinate() == Some(1)`. CDS-anchored intronic positions keep the same
/// primary coordinate plus a signed secondary offset, e.g. `c.-106+2` becomes
/// `coordinate() == Some(-106)` and `offset() == Some(2)`.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{CodingDnaOutcome, NucleotideAnchor, VariantDescription, parse_hgvs};
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let five_prime = parse_hgvs("NM_007373.4:c.-1C>T")?;
/// let three_prime = parse_hgvs("NM_001272071.2:c.*1C>T")?;
///
/// let VariantDescription::CodingDna(CodingDnaOutcome::Known(five_prime_edit)) =
///     five_prime.description
/// else {
///     panic!("expected a coding-DNA edit");
/// };
///
/// let VariantDescription::CodingDna(CodingDnaOutcome::Known(three_prime_edit)) =
///     three_prime.description
/// else {
///     panic!("expected a coding-DNA edit");
/// };
///
/// assert_eq!(
///     five_prime_edit.location.start().unwrap().anchor(),
///     Some(NucleotideAnchor::RelativeCdsStart)
/// );
/// assert_eq!(five_prime_edit.location.start().unwrap().coordinate(), Some(-1));
///
/// assert_eq!(
///     three_prime_edit.location.start().unwrap().anchor(),
///     Some(NucleotideAnchor::RelativeCdsEnd)
/// );
/// assert_eq!(three_prime_edit.location.start().unwrap().coordinate(), Some(1));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleotideCoordinate {
    /// Known nucleotide coordinate with anchor and optional offset, such as
    /// `357+1`, `-1`, or `*1`.
    Known {
        anchor: NucleotideAnchor,
        coordinate: i32,
        offset: i32,
    },
    /// Unknown nucleotide coordinate written as `?`.
    Unknown,
}

impl NucleotideCoordinate {
    /// Builds a known nucleotide coordinate.
    pub fn known(anchor: NucleotideAnchor, coordinate: i32, offset: i32) -> Self {
        Self::Known {
            anchor,
            coordinate,
            offset,
        }
    }

    /// Returns `true` when this coordinate is known.
    pub fn is_known(&self) -> bool {
        matches!(self, Self::Known { .. })
    }

    /// Returns `true` when this coordinate is written as `?`.
    pub fn is_unknown(&self) -> bool {
        matches!(self, Self::Unknown)
    }

    /// Returns the anchor when this coordinate is known.
    pub fn anchor(&self) -> Option<NucleotideAnchor> {
        match self {
            Self::Known { anchor, .. } => Some(*anchor),
            Self::Unknown => None,
        }
    }

    /// Returns the primary coordinate when it is known.
    pub fn coordinate(&self) -> Option<i32> {
        match self {
            Self::Known { coordinate, .. } => Some(*coordinate),
            Self::Unknown => None,
        }
    }

    /// Returns the offset when this coordinate is known.
    pub fn offset(&self) -> Option<i32> {
        match self {
            Self::Known { offset, .. } => Some(*offset),
            Self::Unknown => None,
        }
    }

    /// Returns `true` for intronic coordinates such as `357+1`, `-106+2`,
    /// and `*639-1`.
    pub fn is_intronic(&self) -> bool {
        self.offset().map_or(false, |offset| offset != 0)
    }

    /// Returns `true` if variant's location is relative to the CDS start, such
    /// as `c.-1` and `c.-106+2`.
    pub fn is_cds_start_anchored(&self) -> bool {
        matches!(self.anchor(), Some(NucleotideAnchor::RelativeCdsStart))
    }

    /// Returns `true` if variant's location is relative to the CDS end, such
    /// as `c.*1` and `c.*639-1`.
    pub fn is_cds_end_anchored(&self) -> bool {
        matches!(self.anchor(), Some(NucleotideAnchor::RelativeCdsEnd))
    }

    /// Returns `true` for exonic 5' UTR coordinates such as `c.-81`.
    pub fn is_five_prime_utr(&self) -> bool {
        self.is_cds_start_anchored() && self.offset() == Some(0)
    }

    /// Returns `true` for exonic 3' UTR coordinates such as `c.*24`.
    pub fn is_three_prime_utr(&self) -> bool {
        self.is_cds_end_anchored() && self.offset() == Some(0)
    }
}
