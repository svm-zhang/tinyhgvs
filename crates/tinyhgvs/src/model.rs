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
///         assert_eq!(description.location.start().unwrap().coordinate, -1);
///     }
///     _ => unreachable!("expected nucleotide variant"),
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
    /// DNA or RNA allele container with one written allele, an optional second
    /// established allele, and any later unphased additions.
    NucleotideAllele(AlleleVariant<NucleotideVariant>),
    /// Protein allele container with one written allele, an optional second
    /// established allele, and any later unphased additions.
    ProteinAllele(AlleleVariant<ProteinVariant>),
    Protein(ProteinVariant),
}

/// Phase relationship between two established alleles.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{AllelePhase, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NC_000001.11:g.123G>A(;)345del").unwrap();
///
/// match variant.description {
///     VariantDescription::NucleotideAllele(allele) => {
///         assert_eq!(allele.phase, Some(AllelePhase::Uncertain));
///     }
///     _ => unreachable!("expected nucleotide allele"),
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AllelePhase {
    Trans,
    Uncertain,
}

/// One allele containing one or more inner variants.
///
/// Variants inside one allele are implicitly written in cis.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NC_000001.11:g.[123G>A;345del]").unwrap();
///
/// match variant.description {
///     VariantDescription::NucleotideAllele(allele) => {
///         assert_eq!(allele.allele_one.variants.len(), 2);
///     }
///     _ => unreachable!("expected nucleotide allele"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele<T> {
    pub variants: Vec<T>,
}

impl<T> Allele<T> {
    /// Builds one allele from all its carrying variants.
    pub fn from_variants(variants: Vec<T>) -> Self {
        Self { variants }
    }

    /// Returns the inner variants carried by this allele.
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.variants.iter()
    }
}

impl<'a, T> IntoIterator for &'a Allele<T> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.variants.iter()
    }
}

/// Allele container holding an initial allele, an optional second established
/// allele, and any later unphased alleles.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{AllelePhase, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_004006.2:c.[2376G>C];[2376=]").unwrap();
///
/// match variant.description {
///     VariantDescription::NucleotideAllele(allele) => {
///         assert_eq!(allele.allele_one.variants.len(), 1);
///         assert!(allele.allele_two.is_some());
///         assert_eq!(allele.phase, Some(AllelePhase::Trans));
///         assert_eq!(allele.iter().count(), 2);
///     }
///     _ => unreachable!("expected nucleotide allele"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleVariant<T> {
    pub allele_one: Allele<T>,
    pub allele_two: Option<Allele<T>>,
    pub phase: Option<AllelePhase>,
    pub alleles_unphased: Vec<Allele<T>>,
}

impl<T> AlleleVariant<T> {
    /// Returns all written alleles in order.
    pub fn iter(&self) -> impl Iterator<Item = &Allele<T>> {
        std::iter::once(&self.allele_one)
            .chain(self.allele_two.iter())
            .chain(self.alleles_unphased.iter())
    }

    /// Returns the established trans allele pair when the relation is known.
    pub fn phased_alleles(&self) -> Option<(&Allele<T>, &Allele<T>)> {
        match (self.phase, self.allele_two.as_ref()) {
            (Some(AllelePhase::Trans), Some(allele_two)) => Some((&self.allele_one, allele_two)),
            _ => None,
        }
    }

    /// Returns any later alleles written in uncertain relation to the
    /// established allele state.
    pub fn unphased_alleles(&self) -> &[Allele<T>] {
        &self.alleles_unphased
    }
}

impl<'a, T> IntoIterator for &'a AlleleVariant<T> {
    type Item = &'a Allele<T>;
    type IntoIter = std::iter::Chain<
        std::iter::Chain<std::iter::Once<&'a Allele<T>>, std::option::Iter<'a, Allele<T>>>,
        std::slice::Iter<'a, Allele<T>>,
    >;

    fn into_iter(self) -> Self::IntoIter {
        std::iter::once(&self.allele_one)
            .chain(self.allele_two.iter())
            .chain(self.alleles_unphased.iter())
    }
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
///         assert_eq!(description.location.start().unwrap().coordinate, 357);
///         assert_eq!(description.location.start().unwrap().offset, 1);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideVariant {
    pub location: Location<NucleotideCoordinate>,
    pub edit: NucleotideEdit,
}

/// Parsed protein consequence.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.1:p.(Trp24Ter)").unwrap();
/// let extension = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(description) => {
///         assert!(description.is_predicted);
///         assert!(matches!(description.effect, ProteinEffect::Edit { .. }));
///     }
///     _ => unreachable!("expected protein variant"),
/// }
///
/// match extension.description {
///     VariantDescription::Protein(description) => match description.effect {
///         ProteinEffect::Edit { edit: ProteinEdit::Extension(_), .. } => {}
///         _ => unreachable!("expected protein extension"),
///     },
///     _ => unreachable!("expected protein variant"),
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
        location: Location<ProteinCoordinate>,
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

/// Protein terminus toward which an extension variant extends.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, ProteinExtensionTerminal, VariantDescription, parse_hgvs};
///
/// let n_terminal = parse_hgvs("NP_003997.2:p.Met1ext-5").unwrap();
/// let c_terminal = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17").unwrap();
///
/// let extract_terminal = |variant: tinyhgvs::HgvsVariant| match variant.description {
///     VariantDescription::Protein(description) => match description.effect {
///         ProteinEffect::Edit { edit: ProteinEdit::Extension(extension), .. } => {
///             extension.to_terminal
///         }
///         _ => unreachable!("expected protein extension"),
///     },
///     _ => unreachable!("expected protein variant"),
/// };
///
/// assert_eq!(extract_terminal(n_terminal), ProteinExtensionTerminal::N);
/// assert_eq!(extract_terminal(c_terminal), ProteinExtensionTerminal::C);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinExtensionTerminal {
    N,
    C,
}

/// Model describing a protein extension consequence.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, VariantDescription, parse_hgvs};
///
/// let n_terminal = parse_hgvs("NP_003997.2:p.Met1ext-5").unwrap();
/// let c_terminal = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17").unwrap();
///
/// let extract_extension = |variant: tinyhgvs::HgvsVariant| match variant.description {
///     VariantDescription::Protein(description) => match description.effect {
///         ProteinEffect::Edit { edit: ProteinEdit::Extension(extension), .. } => extension,
///         _ => unreachable!("expected protein extension"),
///     },
///     _ => unreachable!("expected protein variant"),
/// };
///
/// let n_terminal = extract_extension(n_terminal);
/// assert!(n_terminal.to_residue.is_none());
/// assert_eq!(n_terminal.terminal_ordinal, Some(-5));
///
/// let c_terminal = extract_extension(c_terminal);
/// assert_eq!(c_terminal.to_residue.as_deref(), Some("Gln"));
/// assert_eq!(c_terminal.terminal_ordinal, Some(17));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinExtensionEdit {
    pub to_terminal: ProteinExtensionTerminal,
    pub to_residue: Option<String>,
    pub terminal_ordinal: Option<i32>,
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
///         assert_eq!(location.start().unwrap().residue, "Lys");
///         assert_eq!(location.end().unwrap().residue, "Val");
///     }
///     _ => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval<T> {
    pub start: T,
    pub end: Option<T>,
}

/// Main edited location on a nucleotide or protein variant/effect.
///
/// Known locations keep the current one-level interval shape. Uncertain
/// locations wrap the left and right uncertain regions as intervals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Location<T> {
    Known(Interval<T>),
    Uncertain(Interval<Interval<T>>),
}

impl<T> Location<T> {
    /// Builds one known location from a plain interval.
    pub fn from_known(value: Interval<T>) -> Self {
        Self::Known(value)
    }

    /// Builds one uncertain location from uncertain left/right regions.
    pub fn from_uncertain(value: Interval<Interval<T>>) -> Self {
        Self::Uncertain(value)
    }

    /// Returns `true` when the location is written with uncertain-region
    /// syntax.
    pub fn is_uncertain(&self) -> bool {
        matches!(self, Self::Uncertain(_))
    }

    /// Returns `true` for one known position.
    pub fn is_pos(&self) -> bool {
        matches!(self, Self::Known(interval) if interval.end.is_none())
    }

    /// Returns `true` for interval-shaped locations, whether known or
    /// uncertain.
    pub fn is_interval(&self) -> bool {
        !self.is_pos()
    }

    /// Returns the left known position when the location is known.
    pub fn start(&self) -> Option<&T> {
        match self {
            Self::Known(interval) => Some(&interval.start),
            Self::Uncertain(_) => None,
        }
    }

    /// Returns the right known position when the location is a known interval.
    pub fn end(&self) -> Option<&T> {
        match self {
            Self::Known(interval) => interval.end.as_ref(),
            Self::Uncertain(_) => None,
        }
    }

    /// Returns the left uncertain region when the location is uncertain.
    pub fn l_interval(&self) -> Option<&Interval<T>> {
        match self {
            Self::Known(_) => None,
            Self::Uncertain(interval) => Some(&interval.start),
        }
    }

    /// Returns the right uncertain region when the location is an uncertain
    /// interval.
    pub fn r_interval(&self) -> Option<&Interval<T>> {
        match self {
            Self::Known(_) => None,
            Self::Uncertain(interval) => interval.end.as_ref(),
        }
    }
}

/// Primary nucleotide coordinate value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateKind {
    Known(i32),
    Unknown,
}

impl CoordinateKind {
    /// Returns the numeric coordinate when it is known.
    pub fn value(self) -> Option<i32> {
        match self {
            Self::Known(value) => Some(value),
            Self::Unknown => None,
        }
    }
}

impl PartialEq<i32> for CoordinateKind {
    fn eq(&self, other: &i32) -> bool {
        matches!(self, Self::Known(value) if value == other)
    }
}

impl PartialEq<CoordinateKind> for i32 {
    fn eq(&self, other: &CoordinateKind) -> bool {
        other == self
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

/// Nucleotide coordinate with explicit anchor and offset.
///
/// The primary coordinate keeps the sign written in the HGVS string. For
/// example, `c.-1` becomes `coordinate == -1`, while `c.*1` becomes
/// `coordinate == 1`. CDS-anchored intronic positions keep the same primary
/// coordinate plus a signed secondary offset, e.g. `c.-106+2` becomes
/// `coordinate == -106` and `offset == 2`.
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
///         assert_eq!(description.location.start().unwrap().anchor, NucleotideAnchor::RelativeCdsStart);
///         assert_eq!(description.location.start().unwrap().coordinate, -1);
///         assert_eq!(description.location.start().unwrap().offset, 0);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
///
/// match three_prime.description {
///     VariantDescription::Nucleotide(description) => {
///         assert_eq!(description.location.start().unwrap().anchor, NucleotideAnchor::RelativeCdsEnd);
///         assert_eq!(description.location.start().unwrap().coordinate, 1);
///         assert_eq!(description.location.start().unwrap().offset, 0);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
///
/// match five_prime_intronic.description {
///     VariantDescription::Nucleotide(description) => {
///         assert_eq!(description.location.start().unwrap().anchor, NucleotideAnchor::RelativeCdsStart);
///         assert_eq!(description.location.start().unwrap().coordinate, -106);
///         assert_eq!(description.location.start().unwrap().offset, 2);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideCoordinate {
    /// Anchor type used to describe `coordinate`.
    pub anchor: NucleotideAnchor,
    /// Primary coordinate written in the HGVS string.
    pub coordinate: CoordinateKind,
    /// Secondary displacement written after the primary coordinate.
    pub offset: i32,
}

impl NucleotideCoordinate {
    /// Returns `true` when the primary coordinate is written as `?`.
    pub fn is_unknown(&self) -> bool {
        matches!(self.coordinate, CoordinateKind::Unknown)
    }

    /// Returns `true` for intronic coordinates such as `357+1`, `-106+2`,
    /// and `*639-1`.
    pub fn is_intronic(&self) -> bool {
        self.offset != 0
    }

    /// Returns `true` if variant's location is relative to the CDS start, such
    /// as `c.-1` and `c.-106+2`.
    pub fn is_cds_start_anchored(&self) -> bool {
        matches!(self.anchor, NucleotideAnchor::RelativeCdsStart)
    }

    /// Returns `true` if variant's location is relative to the CDS end, such
    /// as `c.*1` and `c.*639-1`.
    pub fn is_cds_end_anchored(&self) -> bool {
        matches!(self.anchor, NucleotideAnchor::RelativeCdsEnd)
    }

    /// Returns `true` for exonic 5' UTR coordinates such as `c.-81`.
    pub fn is_five_prime_utr(&self) -> bool {
        self.is_cds_start_anchored() && self.offset == 0
    }

    /// Returns `true` for exonic 3' UTR coordinates such as `c.*24`.
    pub fn is_three_prime_utr(&self) -> bool {
        self.is_cds_end_anchored() && self.offset == 0
    }
}

impl Interval<NucleotideCoordinate> {
    fn is_end_bound_unknown(&self) -> bool {
        self.end
            .as_ref()
            .map_or(false, NucleotideCoordinate::is_unknown)
    }

    /// Returns `true` when either bound of the interval is unknown, i.e. `?_B`,
    /// `A_?`, `?_?`
    pub fn has_unknown_bound(&self) -> bool {
        self.start.is_unknown() || self.is_end_bound_unknown()
    }

    /// Returns `true` when both sides of the interval are unknown, i.e. ?_?
    pub fn is_fully_unknown(&self) -> bool {
        self.start.is_unknown() && self.is_end_bound_unknown()
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
///         assert_eq!(location.start().unwrap().residue, "Trp");
///         assert_eq!(location.start().unwrap().ordinal, 24);
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
///     _ => unreachable!("expected nucleotide variant"),
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
    /// Protein extension such as `p.Met1ext-5` or `p.Ter110GlnextTer17`.
    Extension(ProteinExtensionEdit),
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
///     _ => unreachable!("expected protein variant"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinSequence {
    pub residues: Vec<String>,
}
