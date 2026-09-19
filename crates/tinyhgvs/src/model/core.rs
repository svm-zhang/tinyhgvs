use super::{AlleleForm, NucleotideEdit, ProteinEdit};

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
///         assert_eq!(description.location.start().unwrap().coordinate().unwrap(), -1);
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

impl HgvsVariant {
    pub fn from(
        reference: Option<ReferenceSpec>,
        coordinate_system: CoordinateSystem,
        description: VariantDescription,
    ) -> Self {
        Self {
            reference,
            coordinate_system,
            description,
        }
    }
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
    CodingDna,
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
            Self::CodingDna => "c",
            Self::Rna => "r",
            Self::Protein => "p",
        }
    }
}

/// Top-level variant description for nucleotide or protein syntax.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantDescription {
    Genomic(GenomicOutcome),
    Rna(RnaOutcome),
    CodingDna(CodingDnaOutcome),
    Protein(ProteinOutcome),

    GenomicAllele(AlleleForm<GenomicOutcome>),
    CodingDnaAllele(AlleleForm<CodingDnaOutcome>),
    RnaAllele(AlleleForm<RnaOutcome>),
    ProteinAllele(AlleleForm<ProteinOutcome>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GenomicOutcome {
    Known(NucleotideEdit),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CodingDnaOutcome {
    Known(NucleotideEdit),
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RnaOutcome {
    // r.A, r.(A)
    Produced {
        edit: NucleotideEdit,
        certainty: OutcomeCertainty,
    },
    // =, (=)
    NoChange(OutcomeCertainty),
    // r.0, r.0?
    NoneProduced(OutcomeCertainty),
    // r.spl, r.spl?
    UncertainSplicing,
    // r.?
    Unknown,
    // r.(?) - Known to exist, but content is unknown
    Indeterminate,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutcomeCertainty {
    Certain,
    Predicted,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinOutcome {
    // p.?
    Unknown,
    // p.0, p.0?
    NoneProduced(OutcomeCertainty),
    Produced {
        edit: ProteinEdit,
        certainty: OutcomeCertainty,
    },
}

impl From<NucleotideEdit> for GenomicOutcome {
    fn from(edit: NucleotideEdit) -> Self {
        Self::Known(edit)
    }
}

impl From<NucleotideEdit> for CodingDnaOutcome {
    fn from(edit: NucleotideEdit) -> Self {
        Self::Known(edit)
    }
}

impl From<NucleotideEdit> for RnaOutcome {
    fn from(edit: NucleotideEdit) -> Self {
        Self::Produced {
            edit,
            certainty: OutcomeCertainty::Certain,
        }
    }
}

impl From<ProteinEdit> for ProteinOutcome {
    fn from(edit: ProteinEdit) -> Self {
        Self::Produced {
            edit,
            certainty: OutcomeCertainty::Certain,
        }
    }
}
