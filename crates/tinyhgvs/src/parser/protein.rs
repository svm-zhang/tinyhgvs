//! Protein outcome, edit, and allele parsers.

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::char;
use nom::combinator::{map, map_res, value};
use nom::multi::{many0, many1, separated_list1};
use nom::sequence::{delimited, pair, preceded, separated_pair};
use nom::Parser;

use super::core::{parse_i32, parse_quantity, range_with};
use super::repeat::{known_repeat_copy, uncertain_repeat_copy};
use super::ParseResult;
use crate::model::{
    Allele, AlleleForm, AllelePhase, AlleleVariant, DerivedAllele, Interval, Location,
    OutcomeCertainty, ProteinCoordinate, ProteinEdit, ProteinEditKind, ProteinExtensionEdit,
    ProteinExtensionTerminal, ProteinFrameshiftStop, ProteinFrameshiftStopKind, ProteinOutcome,
    ProteinSequence, RepeatEdit, VariantDescription,
};

const PROTEIN_SYMBOLS: &[&str] = &[
    "Ter", "Sec", "Pyl", "Xaa", "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His",
    "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "*", "A", "R",
    "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
];

/// Parses a produced protein outcome.
///
/// Examples: `Trp24Ter`, `(Trp24Ter)`
pub(super) fn protein_outcome(input: &str) -> ParseResult<'_, ProteinOutcome> {
    alt((
        // (Ser68Arg)
        map(delimited(char('('), protein_edit, char(')')), |edit| {
            ProteinOutcome::Produced {
                edit,
                certainty: OutcomeCertainty::Predicted,
            }
        }),
        // Ser68Arg
        map(protein_edit, |edit| ProteinOutcome::Produced {
            edit,
            certainty: OutcomeCertainty::Certain,
        }),
    ))
    .parse(input)
}

/// Parses special protein outcomes.
///
/// Examples: `?`, `0`, `0?`
pub(super) fn special_protein_outcome(input: &str) -> ParseResult<'_, ProteinOutcome> {
    alt((
        // p.?
        value(ProteinOutcome::Unknown, char('?')),
        // p.0?
        value(
            ProteinOutcome::NoneProduced(OutcomeCertainty::Predicted),
            tag("0?"),
        ),
        // p.0
        value(
            ProteinOutcome::NoneProduced(OutcomeCertainty::Certain),
            char('0'),
        ),
    ))
    .parse(input)
}

/// Parses one or more protein outcomes written on the same allele.
///
/// Examples: `Ser68Arg;Asn594del`, `(Ser68Arg;Asn594del)`
pub(super) fn protein_variants_on_allele(input: &str) -> ParseResult<'_, Vec<ProteinOutcome>> {
    alt((
        // (Ser68Arg;Asn594del)
        map(
            delimited(
                char('('),
                separated_list1(char(';'), protein_edit),
                char(')'),
            ),
            |edits| {
                edits
                    .into_iter()
                    .map(|edit| ProteinOutcome::Produced {
                        edit,
                        certainty: OutcomeCertainty::Predicted,
                    })
                    .collect()
            },
        ),
        // Ser68Arg;Asn594del
        // Phe233Leu;(Cys690Trp)
        separated_list1(char(';'), protein_outcome),
    ))
    .parse(input)
}

/// Parses one protein allele component.
///
/// Examples: `[?]`, `[0]`, `[Ser68Arg]`, `[Ser68Arg;Asn594del]`
pub(super) fn protein_allele_component(input: &str) -> ParseResult<'_, Allele<ProteinOutcome>> {
    map(
        delimited(
            char('['),
            alt((
                // [?]
                map(char('?'), |_| vec![ProteinOutcome::Unknown]),
                // [0]
                map(char('0'), |_| {
                    vec![ProteinOutcome::NoneProduced(OutcomeCertainty::Certain)]
                }),
                protein_variants_on_allele,
            )),
            char(']'),
        ),
        Allele::from_variants,
    )
    .parse(input)
}

/// Parses one in-cis protein allele.
///
/// Example: `[Ser68Arg;Asn594del]`
pub(super) fn protein_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    map(protein_allele_component, |allele| AlleleVariant {
        allele_one: allele,
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

/// Parses protein alleles in trans.
///
/// Example: `[Ser68Arg];[Ser68=]`
pub(super) fn protein_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    map(
        separated_pair(
            protein_allele_component,
            char(';'),
            protein_allele_component,
        ),
        |(a1, a2)| AlleleVariant {
            allele_one: a1,
            allele_two: Some(a2),
            phase: Some(AllelePhase::Trans),
            variants_unphased: vec![],
        },
    )
    .parse(input)
}

/// Parses protein alleles with uncertain phase.
///
/// Example: `Ser68Arg(;)Asn594del`
pub(super) fn protein_uncertain_allele(
    input: &str,
) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    map(
        separated_pair(protein_outcome, tag("(;)"), protein_outcome),
        |(a1, a2)| AlleleVariant {
            allele_one: Allele::from_variants(vec![a1]),
            allele_two: Some(Allele::from_variants(vec![a2])),
            phase: Some(AllelePhase::Uncertain),
            variants_unphased: vec![],
        },
    )
    .parse(input)
}

/// Parses supported single protein allele forms.
///
/// This is the main entrance for handling possible protein allele descriptions.
///
/// Examples: `[A;B]`, `[A];[B]`, `A(;)B`
pub(super) fn protein_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    alt((
        // [Ser68Arg];[Ser68=]
        protein_trans_allele,
        // Ser68Arg(;)Asn594del
        protein_uncertain_allele,
        // [Ser68Arg;Asn594del]
        protein_cis_allele,
    ))
    .parse(input)
}

/// Parses a derived protein allele form.
///
/// Example: `[Lys31Asn,Val25_Lys31del,Ser68Arg]`
pub(super) fn protein_derived_allele_form(
    input: &str,
) -> ParseResult<'_, DerivedAllele<ProteinOutcome>> {
    // [Lys31Asn,Val25_Lys31del,Ser68Arg]
    let (input, _) = char('[')(input)?;

    let (input, first) = protein_outcome(input)?;
    let (input, _) = char(',')(input)?;
    let (input, second) = protein_outcome(input)?;
    let (input, rest) = many0(preceded(char(','), protein_outcome)).parse(input)?;

    let (input, _) = char(']')(input)?;

    let mut outcomes = Vec::with_capacity(2 + rest.len());
    outcomes.push(first);
    outcomes.push(second);
    outcomes.extend(rest);

    Ok((input, DerivedAllele::from_outcomes(outcomes)))
}

/// Parses alternative protein allele forms.
///
/// Example: `[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]`
pub(super) fn protein_alternative_allele_form(
    input: &str,
) -> ParseResult<'_, Vec<AlleleVariant<ProteinOutcome>>> {
    fn alternate(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
        alt((
            // [(A)(;)(B)]
            delimited(char('['), protein_uncertain_allele, char(']')),
            // [A], [(A)], [A;B], ...
            protein_cis_allele,
        ))
        .parse(input)
    }

    let (input, first) = alternate(input)?;
    let (input, _) = char('^')(input)?;
    let (input, second) = alternate(input)?;

    let (input, rest) = many0(preceded(char('^'), alternate)).parse(input)?;

    let mut alternatives = Vec::with_capacity(2 + rest.len());
    alternatives.push(first);
    alternatives.push(second);
    alternatives.extend(rest);

    Ok((input, alternatives))
}

/// Parses a protein (p.) description.
///
/// Examples: `p.Trp24Ter`, `p.?`, `p.[Ser68Arg];[Ser68=]`
pub(super) fn protein_description(input: &str) -> ParseResult<'_, VariantDescription> {
    alt((
        // p.[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]
        map(protein_alternative_allele_form, |alternatives| {
            VariantDescription::ProteinAllele(AlleleForm::Alternative(alternatives))
        }),
        // p.[Lys31Asn,Val25_Lys31del,Ser68Arg]
        map(protein_derived_allele_form, |derived| {
            VariantDescription::ProteinAllele(AlleleForm::Derived(derived))
        }),
        // p.[Ser68Arg];[Ser68=]
        map(protein_allele, |variant| {
            VariantDescription::ProteinAllele(AlleleForm::Single(variant))
        }),
        // p.?, p.0, p.0?
        map(special_protein_outcome, VariantDescription::Protein),
        // p.Trp24Ter, p.(Trp24Ter)
        map(protein_outcome, VariantDescription::Protein),
    ))
    .parse(input)
}

/// Parses one protein edit: location + edit.
///
/// Examples: `Trp24Ter`, `Arg97ProfsTer23`, `Ter110GlnextTer17`
pub(super) fn protein_edit(input: &str) -> ParseResult<'_, ProteinEdit> {
    map_res(
        pair(protein_location, protein_edit_kind),
        build_protein_edit_effect,
    )
    .parse(input)
}

/// Builds a protein edit after validating edit/location combinations.
///
/// This rejects malformed extension combinations such as C-terminal extension
/// syntax on a non-terminating residue.
pub(super) fn build_protein_edit_effect(
    (location, kind): (Location<ProteinCoordinate>, ProteinEditKind),
) -> Result<ProteinEdit, ()> {
    let location = resolve_protein_effect_location(&location, &kind).ok_or(())?;
    Ok(ProteinEdit { location, kind })
}

/// Resolves and validates the protein location for extension edits.
///
/// Early returns reject extension syntax that cannot apply to the parsed
/// location, such as interval locations or N-terminal extension away from
/// `Met1`.
pub(super) fn resolve_protein_effect_location(
    location: &Location<ProteinCoordinate>,
    edit: &ProteinEditKind,
) -> Option<Location<ProteinCoordinate>> {
    let ProteinEditKind::Extension(extension) = edit else {
        return Some(location.clone());
    };

    let Location::Known(location) = location else {
        return None;
    };

    if location.end.is_some() {
        return None;
    }

    let mut start = location.start.clone();

    match extension.to_terminal {
        ProteinExtensionTerminal::N => {
            if start.residue != "Met"
                || start.ordinal != 1
                || extension.to_residue.is_some()
                || !matches!(extension.terminal_ordinal, Some(ordinal) if ordinal < 0)
            {
                return None;
            }
        }
        ProteinExtensionTerminal::C => {
            if start.residue != "Ter"
                || extension.to_residue.is_none()
                || matches!(extension.terminal_ordinal, Some(ordinal) if ordinal <= 0)
            {
                return None;
            }
            start.residue = "Ter".to_string();
        }
    }

    Some(Location::from_known(Interval { start, end: None }))
}

/// Parses various protein edit families: duplication, deletion, etc.
///
/// Examples: `Ter`, `del`, `dup`, `insAla`, `fs`, `GlnextTer17`
pub(super) fn protein_edit_kind(input: &str) -> ParseResult<'_, ProteinEditKind> {
    alt((
        // (=)
        value(
            ProteinEditKind::NoChange(OutcomeCertainty::Predicted),
            tag("(=)"),
        ),
        // =
        value(
            ProteinEditKind::NoChange(OutcomeCertainty::Certain),
            char('='),
        ),
        // delinsGly
        map(preceded(tag("delins"), protein_sequence), |sequence| {
            ProteinEditKind::DeletionInsertion { sequence }
        }),
        // del
        value(ProteinEditKind::Deletion, tag("del")),
        // dup
        value(ProteinEditKind::Duplication, tag("dup")),
        // [10], [(70_80)]
        protein_repeat,
        // ext-5, GlnextTer17
        protein_extension_edit,
        // fs, ProfsTer23
        protein_frameshift_edit,
        // insAla
        map(preceded(tag("ins"), protein_sequence), |sequence| {
            ProteinEditKind::Insertion { sequence }
        }),
        // Ter, Asp
        map(protein_symbol, |to| ProteinEditKind::Substitution { to }),
    ))
    .parse(input)
}

/// Parses one protein location, known or uncertain.
///
/// Examples: `Trp24`, `Lys23_Val25`, `(Ala123_Pro131)`
pub(super) fn protein_location(input: &str) -> ParseResult<'_, Location<ProteinCoordinate>> {
    alt((
        // (Ala123_Pro131) and (Ala123_Pro131)_(Gly140_Leu142)
        map(protein_uncertain_location, Location::from_uncertain),
        // Trp24 and Lys23_Val25
        map(protein_interval, Location::from_known),
    ))
    .parse(input)
}

/// Parses a single protein position or a known interval.
///
/// Examples: `Ala237`, `Ala237_Pro161`
pub(super) fn protein_interval(input: &str) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    alt((
        |input| range_with(input, protein_coordinate),
        map(protein_coordinate, |start| Interval { start, end: None }),
    ))
    .parse(input)
}

/// Parses one protein uncertain interval unit.
///
/// Example: `(Ala237_Pro161)`
pub(super) fn protein_uncertain_interval(
    input: &str,
) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    delimited(char('('), protein_interval, char(')')).parse(input)
}

/// Parses protein locations written with uncertain-region syntax.
///
/// Examples: `(Ala237_Pro161)`, `(Ala237_Pro161)_(Gly170_Leu180)`
pub(super) fn protein_uncertain_location(
    input: &str,
) -> ParseResult<'_, Interval<Interval<ProteinCoordinate>>> {
    alt((
        |input| range_with(input, protein_uncertain_interval),
        map(protein_uncertain_interval, |start| Interval {
            start,
            end: None,
        }),
    ))
    .parse(input)
}

/// Parses a protein symbol followed by its ordinal.
///
/// Examples: `Trp24`, `Ala237`
pub(super) fn protein_coordinate(input: &str) -> ParseResult<'_, ProteinCoordinate> {
    map(pair(protein_symbol, parse_i32), |(residue, ordinal)| {
        ProteinCoordinate { residue, ordinal }
    })
    .parse(input)
}

/// Parses protein repeat copy syntax after a protein location.
///
/// Examples: `[10]`, `[(70_80)]`
pub(super) fn protein_repeat(input: &str) -> ParseResult<'_, ProteinEditKind> {
    // p.Ala2[10]
    // p.(Gln18)[(70_80)]
    map(
        alt((known_repeat_copy, uncertain_repeat_copy)),
        |quantity| {
            ProteinEditKind::Repeat(RepeatEdit {
                unit: None,
                quantity,
            })
        },
    )
    .parse(input)
}

/// Parses N-terminal and C-terminal protein extension syntax.
///
/// Examples: `ext-5`, `GlnextTer17`
pub(super) fn protein_extension_edit(input: &str) -> ParseResult<'_, ProteinEditKind> {
    alt((
        map(
            preceded(tag("ext"), protein_n_terminal_extension_ordinal),
            |terminal_ordinal| {
                ProteinEditKind::Extension(ProteinExtensionEdit {
                    to_terminal: ProteinExtensionTerminal::N,
                    to_residue: None,
                    terminal_ordinal: Some(terminal_ordinal),
                })
            },
        ),
        map(
            pair(
                protein_extension_residue,
                protein_c_terminal_extension_state,
            ),
            |(to_residue, terminal_ordinal)| {
                ProteinEditKind::Extension(ProteinExtensionEdit {
                    to_terminal: ProteinExtensionTerminal::C,
                    to_residue: Some(to_residue),
                    terminal_ordinal,
                })
            },
        ),
    ))
    .parse(input)
}

/// Parses the required negative ordinal in N-terminal extension syntax.
///
/// Example: `-5`
pub(super) fn protein_n_terminal_extension_ordinal(input: &str) -> ParseResult<'_, i32> {
    map(preceded(char('-'), parse_i32), |ordinal| -ordinal).parse(input)
}

/// Parses the residue replacing the reference stop codon in C-terminal extension syntax.
///
/// Example: `Gln`
pub(super) fn protein_extension_residue(input: &str) -> ParseResult<'_, String> {
    let (input, residue) = protein_symbol(input)?;

    if residue == "Ter" {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )))
    } else {
        Ok((input, residue))
    }
}

/// Parses the terminal state in C-terminal extension syntax.
///
/// Examples: `extTer17`, `extTer?`
pub(super) fn protein_c_terminal_extension_state(input: &str) -> ParseResult<'_, Option<i32>> {
    preceded(
        tag("ext"),
        alt((
            value(None, pair(alt((tag("Ter"), tag("*"))), char('?'))),
            map(preceded(alt((tag("Ter"), tag("*"))), parse_i32), Some),
        )),
    )
    .parse(input)
}

/// Parses short and long protein frameshift syntax.
///
/// Examples: `fs`, `ProfsTer23`, `ProfsTer?`
pub(super) fn protein_frameshift_edit(input: &str) -> ParseResult<'_, ProteinEditKind> {
    alt((
        map(
            pair(
                protein_frameshift_residue,
                pair(tag("fs"), protein_frameshift_stop),
            ),
            |(to_residue, (_, stop))| ProteinEditKind::Frameshift {
                to_residue: Some(to_residue),
                stop,
            },
        ),
        value(
            ProteinEditKind::Frameshift {
                to_residue: None,
                stop: ProteinFrameshiftStop {
                    ordinal: None,
                    kind: ProteinFrameshiftStopKind::Omitted,
                },
            },
            tag("fs"),
        ),
    ))
    .parse(input)
}

/// Parses the explicit stop-state in long protein frameshift notation.
///
/// Examples: `Ter23`, `Ter?`
pub(super) fn protein_frameshift_stop(input: &str) -> ParseResult<'_, ProteinFrameshiftStop> {
    alt((
        value(
            ProteinFrameshiftStop {
                ordinal: None,
                kind: ProteinFrameshiftStopKind::Unknown,
            },
            pair(alt((tag("Ter"), tag("*"))), char('?')),
        ),
        map(
            preceded(alt((tag("Ter"), tag("*"))), parse_quantity),
            |ordinal| ProteinFrameshiftStop {
                ordinal: Some(ordinal),
                kind: ProteinFrameshiftStopKind::Known,
            },
        ),
    ))
    .parse(input)
}

/// Parses the explicitly written first residue in long protein frameshift syntax.
///
/// Example: `Pro`
pub(super) fn protein_frameshift_residue(input: &str) -> ParseResult<'_, String> {
    let (input, residue) = protein_symbol(input)?;

    if residue == "Ter" {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )))
    } else {
        Ok((input, residue))
    }
}

/// Parses a contiguous protein sequence.
///
/// Examples: `Ala`, `GlnSerLys`
pub(super) fn protein_sequence(input: &str) -> ParseResult<'_, ProteinSequence> {
    map(many1(protein_symbol), |residues| ProteinSequence {
        residues,
    })
    .parse(input)
}

/// Parses one supported amino-acid symbol.
///
/// Examples: `Trp`, `W`, `*`
pub(super) fn protein_symbol(input: &str) -> ParseResult<'_, String> {
    for symbol in PROTEIN_SYMBOLS {
        if let Some(rest) = input.strip_prefix(symbol) {
            return Ok((rest, normalize_protein_symbol(symbol)));
        }
    }

    Err(nom::Err::Error(nom::error::Error::new(
        input,
        nom::error::ErrorKind::Tag,
    )))
}

pub(super) fn normalize_protein_symbol(symbol: &str) -> String {
    if symbol == "*" {
        "Ter".to_string()
    } else {
        symbol.to_string()
    }
}
