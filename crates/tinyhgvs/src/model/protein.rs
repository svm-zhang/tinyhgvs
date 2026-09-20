//! Protein coordinates, edits, and edit-specific model details.

use super::{Location, OutcomeCertainty, RepeatEdit};

/// A protein edit applied at a protein location.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinEdit {
    pub location: Location<ProteinCoordinate>,
    pub kind: ProteinEditKind,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinEditKind {
    // p.=, p.(=)
    NoChange(OutcomeCertainty),
    // p.Trp24Ter
    Substitution {
        to: ResidueChange,
    },
    // p.Lys23_Val25del
    Deletion,
    // p.Ser68_Arg70dup
    Duplication,
    // p.Ala2[10], p.Arg65_Ser67[12]
    Repeat(RepeatEdit),
    // p.Met1ext-5, p.Ter110GlnextTer17
    Extension(ProteinExtensionEdit),
    // p.Arg97fs, p.Arg97ProfsTer23
    Frameshift {
        to_residue: Option<ResidueChange>,
        stop: ProteinFrameshiftStop,
    },
    // p.Val582_Asn583insAla
    Insertion {
        sequence: ProteinInsertionSequence,
    },
    // p.Ser68_Arg70delinsGly
    DeletionInsertion {
        sequence: ProteinSequence,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResidueChange {
    // p.Trp24Ter, p.Arg97ProfsTer23
    Known(String),
    // p.(Gly719Ala^Ser), p.Gly719(Ala^Ser)fsTer23
    Alternative(Vec<String>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinInsertionSequence {
    // p.Val582_Asn583insAla
    Known(ProteinSequence),
    // p.Ser332_Ser333insXaa, p.Arg78_Gly79insXaa[23]
    Unknown { count: usize },
    // p.Gln746_Lys747ins*63
    Terminating { ordinal: usize },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinEditForm {
    // p.Trp24Ter
    Single(ProteinEdit),
    // p.(Gly23GlufsTer7^Gly23CysfsTer26)
    Alternative(Vec<ProteinEdit>),
}

/// Model describing a stop codon is known (long-form), or omitted (short-form),
/// or unknown (not encountered) due to a frameshift event.
///
/// - "Known" or long-form: `p.Arg97ProfsTer23`
/// - "Omitted" or short-form: `p.Arg97fs`
/// - "Unknown" or "not encountered": `p.Arg97ProfsTer?`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinFrameshiftStopKind {
    // p.Arg97fs
    Omitted,
    // p.Arg97ProfsTer?
    Unknown,
    // p.Arg97ProfsTer23
    Known,
}

/// Protein terminus toward which an extension variant extends.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{
///     ProteinEditForm, ProteinEditKind, ProteinExtensionTerminal, ProteinOutcome,
///     VariantDescription, parse_hgvs,
/// };
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let n_terminal = parse_hgvs("NP_003997.2:p.Met1ext-5")?;
/// let c_terminal = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")?;
///
/// let extract_terminal = |variant: tinyhgvs::HgvsVariant| match variant.description {
///     VariantDescription::Protein(ProteinOutcome::Produced {
///         edit: ProteinEditForm::Single(edit),
///         ..
///     }) => match edit.kind {
///         ProteinEditKind::Extension(extension) => extension.to_terminal,
///         _ => unreachable!("expected protein extension"),
///     },
///     _ => unreachable!("expected protein variant"),
/// };
///
/// assert_eq!(extract_terminal(n_terminal), ProteinExtensionTerminal::N);
/// assert_eq!(extract_terminal(c_terminal), ProteinExtensionTerminal::C);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinExtensionTerminal {
    // p.Met1ext-5
    N,
    // p.Ter110GlnextTer17
    C,
}

/// Model describing a protein extension consequence.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{
///     ProteinEditForm, ProteinEditKind, ProteinExtensionTerminal, ProteinOutcome,
///     VariantDescription, parse_hgvs,
/// };
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")?;
///
/// let VariantDescription::Protein(ProteinOutcome::Produced {
///     edit: ProteinEditForm::Single(edit),
///     ..
/// }) = variant.description else {
///     panic!("expected a produced protein outcome");
/// };
///
/// let ProteinEditKind::Extension(extension) = edit.kind else {
///     panic!("expected protein extension");
/// };
///
/// assert_eq!(extension.to_terminal, ProteinExtensionTerminal::C);
/// assert_eq!(extension.to_residue.as_deref(), Some("Gln"));
/// assert_eq!(extension.terminal_ordinal, Some(17));
/// # Ok(())
/// # }
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
/// use tinyhgvs::{
///     ProteinEditForm, ProteinEditKind, ProteinFrameshiftStopKind, ProteinOutcome,
///     VariantDescription, parse_hgvs,
/// };
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")?;
///
/// let VariantDescription::Protein(ProteinOutcome::Produced {
///     edit: ProteinEditForm::Single(edit),
///     ..
/// }) = variant.description else {
///     panic!("expected a produced protein outcome");
/// };
///
/// let ProteinEditKind::Frameshift { stop, .. } = edit.kind else {
///     panic!("expected protein frameshift");
/// };
///
/// assert_eq!(stop.kind, ProteinFrameshiftStopKind::Known);
/// assert_eq!(stop.ordinal, Some(23));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinFrameshiftStop {
    pub ordinal: Option<usize>,
    pub kind: ProteinFrameshiftStopKind,
}

/// Ordered protein insertion or replacement sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinSequence {
    pub residues: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinCoordinate {
    pub residue: String,
    pub ordinal: i32,
}
