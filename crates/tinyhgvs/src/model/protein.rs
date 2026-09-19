use super::{Location, OutcomeCertainty, RepeatEdit};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinEdit {
    pub location: Location<ProteinCoordinate>,
    pub kind: ProteinEditKind,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinEditKind {
    // p.=, p.(=)
    NoChange(OutcomeCertainty),
    Substitution {
        to: String,
    },
    Deletion,
    Duplication,
    /// Top-level repeated sequence such as `p.Ala2[10]` or
    /// `p.Arg65_Ser67[12]`.
    // Repeat {
    //     count: usize,
    // },
    Repeat(RepeatEdit),
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
///         ProteinEffect::Known { edit: ProteinEdit::Extension(extension), .. } => {
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
///         ProteinEffect::Known { edit: ProteinEdit::Extension(extension), .. } => extension,
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
///         ProteinEffect::Known { edit: ProteinEdit::Frameshift { stop, .. }, .. } => stop,
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
///         let ProteinEffect::Known { edit, .. } = description.effect else {
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinCoordinate {
    pub residue: String,
    pub ordinal: i32,
}
