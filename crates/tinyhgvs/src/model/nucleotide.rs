use super::{CoordinateSystem, Interval, Location, ReferenceSpec, RepeatEdit};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideEdit {
    pub location: Location<NucleotideCoordinate>,
    pub kind: NucleotideEditKind,
}

/// Supported nucleotide edit families.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NucleotideEditKind {
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
    /// Top-level repeated sequence such as `g.123CAG[23]`
    Repeat {
        blocks: Vec<RepeatEdit>,
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
    Repeat(RepeatEdit),
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
///         assert_eq!(item.source_location.start.coordinate().unwrap(), 450);
///         assert_eq!(item.source_location.end.as_ref().unwrap().coordinate().unwrap(), 470);
///     }
///     _ => unreachable!("expected nucleotide variant"),
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
/// use tinyhgvs::{NucleotideAnchor, VariantDescription, parse_hgvs};
///
/// let five_prime = parse_hgvs("NM_007373.4:c.-1C>T").unwrap();
/// let three_prime = parse_hgvs("NM_001272071.2:c.*1C>T").unwrap();
/// let five_prime_intronic = parse_hgvs("NM_001385026.1:c.-106+2T>A").unwrap();
///
/// match five_prime.description {
///     VariantDescription::Nucleotide(description) => {
///         assert_eq!(description.location.start().unwrap().anchor(), Some(NucleotideAnchor::RelativeCdsStart));
///         assert_eq!(description.location.start().unwrap().coordinate(), Some(-1));
///         assert_eq!(description.location.start().unwrap().offset(), Some(0));
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
///
/// match three_prime.description {
///     VariantDescription::Nucleotide(description) => {
///         assert_eq!(description.location.start().unwrap().anchor(), Some(NucleotideAnchor::RelativeCdsEnd));
///         assert_eq!(description.location.start().unwrap().coordinate(), Some(1));
///         assert_eq!(description.location.start().unwrap().offset(), Some(0));
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
///
/// match five_prime_intronic.description {
///     VariantDescription::Nucleotide(description) => {
///         assert_eq!(description.location.start().unwrap().anchor(), Some(NucleotideAnchor::RelativeCdsStart));
///         assert_eq!(description.location.start().unwrap().coordinate(), Some(-106));
///         assert_eq!(description.location.start().unwrap().offset(), Some(2));
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleotideCoordinate {
    /// Known nucleotide coordinate with anchor and optional offset.
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
