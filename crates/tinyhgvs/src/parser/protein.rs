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

pub(super) fn protein_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    map(protein_allele_component, |allele| AlleleVariant {
        allele_one: allele,
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

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

pub(super) fn protein_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    alt((
        protein_trans_allele,
        protein_uncertain_allele,
        protein_cis_allele,
    ))
    .parse(input)
}

pub(super) fn protein_derived_allele_form(
    input: &str,
) -> ParseResult<'_, DerivedAllele<ProteinOutcome>> {
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

/// Parser for protein variant and allele description.
pub(super) fn protein_description(input: &str) -> ParseResult<'_, VariantDescription> {
    alt((
        map(protein_alternative_allele_form, |alternatives| {
            VariantDescription::ProteinAllele(AlleleForm::Alternative(alternatives))
        }),
        map(protein_derived_allele_form, |derived| {
            VariantDescription::ProteinAllele(AlleleForm::Derived(derived))
        }),
        map(protein_allele, |variant| {
            VariantDescription::ProteinAllele(AlleleForm::Single(variant))
        }),
        map(special_protein_outcome, VariantDescription::Protein),
        map(protein_outcome, VariantDescription::Protein),
    ))
    .parse(input)
}

pub(super) fn protein_edit(input: &str) -> ParseResult<'_, ProteinEdit> {
    map_res(
        pair(protein_location, protein_edit_kind),
        build_protein_edit_effect,
    )
    .parse(input)
}

pub(super) fn build_protein_edit_effect(
    (location, kind): (Location<ProteinCoordinate>, ProteinEditKind),
) -> Result<ProteinEdit, ()> {
    let location = resolve_protein_effect_location(&location, &kind).ok_or(())?;
    Ok(ProteinEdit { location, kind })
}

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

/// Parses the currently supported protein edit families.
pub(super) fn protein_edit_kind(input: &str) -> ParseResult<'_, ProteinEditKind> {
    alt((
        value(
            ProteinEditKind::NoChange(OutcomeCertainty::Predicted),
            tag("(=)"),
        ),
        value(
            ProteinEditKind::NoChange(OutcomeCertainty::Certain),
            char('='),
        ),
        map(preceded(tag("delins"), protein_sequence), |sequence| {
            ProteinEditKind::DeletionInsertion { sequence }
        }),
        value(ProteinEditKind::Deletion, tag("del")),
        value(ProteinEditKind::Duplication, tag("dup")),
        protein_repeat,
        protein_extension_edit,
        protein_frameshift_edit,
        map(preceded(tag("ins"), protein_sequence), |sequence| {
            ProteinEditKind::Insertion { sequence }
        }),
        map(protein_symbol, |to| ProteinEditKind::Substitution { to }),
    ))
    .parse(input)
}

/// Parses one supported protein location, known or uncertain.
pub(super) fn protein_location(input: &str) -> ParseResult<'_, Location<ProteinCoordinate>> {
    alt((
        // (Ala123_Pro131) and (Ala123_Pro131)_(Gly140_Leu142)
        map(protein_uncertain_location, Location::from_uncertain),
        // Trp24 and Lys23_Val25
        map(protein_interval, Location::from_known),
    ))
    .parse(input)
}

/**
Parses a single protein position or an interval.

- Single position: `Ala237`
- Interval: `Ala237_Pro161`
*/
pub(super) fn protein_interval(input: &str) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    alt((
        |input| range_with(input, protein_coordinate),
        map(protein_coordinate, |start| Interval { start, end: None }),
    ))
    .parse(input)
}

/**
Parse one protein uncertain interval unit (with parenthesis). This is a wrapper
over protein_interval parser.

- `(Ala237_Pro161)`
*/
pub(super) fn protein_uncertain_interval(
    input: &str,
) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    delimited(char('('), protein_interval, char(')')).parse(input)
}

/**
Parses protein locations written with uncertain-region syntax.
*/
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
pub(super) fn protein_coordinate(input: &str) -> ParseResult<'_, ProteinCoordinate> {
    map(pair(protein_symbol, parse_i32), |(residue, ordinal)| {
        ProteinCoordinate { residue, ordinal }
    })
    .parse(input)
}

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
pub(super) fn protein_n_terminal_extension_ordinal(input: &str) -> ParseResult<'_, i32> {
    map(preceded(char('-'), parse_i32), |ordinal| -ordinal).parse(input)
}

/// Parses the residue replacing the reference stop codon in C-terminal extension syntax.
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
pub(super) fn protein_sequence(input: &str) -> ParseResult<'_, ProteinSequence> {
    map(many1(protein_symbol), |residues| ProteinSequence {
        residues,
    })
    .parse(input)
}

/// Parses one supported amino-acid symbol.
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
