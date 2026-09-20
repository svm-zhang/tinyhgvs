//! Protein outcome, edit, and allele parsers.

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{char, digit1};
use nom::combinator::{map, map_res, opt, value};
use nom::multi::{many0, many1, separated_list1};
use nom::sequence::{delimited, pair, preceded, separated_pair};
use nom::Parser;

use super::core::{parse_i32, parse_quantity, range_with};
use super::repeat::{known_repeat_copy, uncertain_repeat_copy};
use super::ParseResult;
use crate::model::{
    Allele, AlleleForm, AllelePhase, AlleleVariant, DerivedAllele, Interval, Location,
    OutcomeCertainty, ProteinCoordinate, ProteinEdit, ProteinEditForm, ProteinEditKind,
    ProteinExtensionEdit, ProteinExtensionTerminal, ProteinFrameshiftStop,
    ProteinFrameshiftStopKind, ProteinInsertionSequence, ProteinOutcome, ProteinSequence,
    RepeatEdit, ResidueChange, VariantDescription,
};

const PROTEIN_SYMBOLS: &[&str] = &[
    "Ter", "Sec", "Pyl", "Xaa", "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His",
    "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "*", "A", "R",
    "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
];

fn protein_alternative_residues(input: &str) -> ParseResult<'_, Vec<String>> {
    let (input, first) = protein_symbol(input)?;
    let (input, _) = char('^')(input)?;
    let (input, second) = protein_symbol(input)?;

    let (input, rest) = many0(preceded(char('^'), protein_symbol)).parse(input)?;

    let mut residues = Vec::with_capacity(2 + rest.len());
    residues.push(first);
    residues.push(second);
    residues.extend(rest);

    Ok((input, residues))
}

fn protein_residue_change(input: &str) -> ParseResult<'_, ResidueChange> {
    alt((
        map(protein_alternative_residues, ResidueChange::Alternative),
        map(protein_symbol, ResidueChange::Known),
    ))
    .parse(input)
}

/// Parses a produced protein outcome.
///
/// Examples: `Trp24Ter`, `(Trp24Ter)`
fn protein_outcome(input: &str) -> ParseResult<'_, ProteinOutcome> {
    alt((
        // (Ser68Arg)
        map(delimited(char('('), protein_edit_form, char(')')), |edit| {
            ProteinOutcome::Produced {
                edit,
                certainty: OutcomeCertainty::Predicted,
            }
        }),
        // Ser68Arg
        map(protein_edit_form, |edit| ProteinOutcome::Produced {
            edit,
            certainty: OutcomeCertainty::Certain,
        }),
    ))
    .parse(input)
}

/// Parses special protein outcomes.
///
/// Examples: `?`, `0`, `0?`
fn special_protein_outcome(input: &str) -> ParseResult<'_, ProteinOutcome> {
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
fn protein_variants_on_allele(input: &str) -> ParseResult<'_, Vec<ProteinOutcome>> {
    alt((
        // (Ser68Arg;Asn594del)
        map(
            delimited(
                char('('),
                separated_list1(char(';'), protein_edit_form),
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
fn protein_allele_component(input: &str) -> ParseResult<'_, Allele<ProteinOutcome>> {
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
fn protein_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
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
fn protein_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
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
fn protein_uncertain_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
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
fn protein_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
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
fn protein_derived_allele_form(input: &str) -> ParseResult<'_, DerivedAllele<ProteinOutcome>> {
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
fn protein_alternative_allele_form(
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
fn protein_edit(input: &str) -> ParseResult<'_, ProteinEdit> {
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
fn build_protein_edit_effect(
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
fn resolve_protein_effect_location(
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
fn protein_edit_kind(input: &str) -> ParseResult<'_, ProteinEditKind> {
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
        // FIXME: insAla
        map(
            preceded(tag("ins"), protein_insertion_sequence),
            |sequence| ProteinEditKind::Insertion { sequence },
        ),
        // FIXME: Ter, Asp
        map(protein_residue_change, |to| ProteinEditKind::Substitution {
            to,
        }),
    ))
    .parse(input)
}

fn protein_alternative_edit_form(input: &str) -> ParseResult<'_, ProteinEditForm> {
    let (input, first) = protein_edit(input)?;
    let (input, _) = char('^')(input)?;
    let (input, second) = protein_edit(input)?;

    let (input, rest) = many0(preceded(char('^'), protein_edit)).parse(input)?;

    let mut edits = Vec::with_capacity(2 + rest.len());
    edits.push(first);
    edits.push(second);
    edits.extend(rest);

    Ok((input, ProteinEditForm::Alternative(edits)))
}

fn protein_edit_form(input: &str) -> ParseResult<'_, ProteinEditForm> {
    alt((
        protein_alternative_edit_form,
        map(protein_edit, ProteinEditForm::Single),
    ))
    .parse(input)
}

/// Parses one protein location, known or uncertain.
///
/// Examples: `Trp24`, `Lys23_Val25`, `(Ala123_Pro131)`
fn protein_location(input: &str) -> ParseResult<'_, Location<ProteinCoordinate>> {
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
fn protein_interval(input: &str) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    alt((
        |input| range_with(input, protein_coordinate),
        map(protein_coordinate, |start| Interval { start, end: None }),
    ))
    .parse(input)
}

/// Parses one protein uncertain interval unit.
///
/// Example: `(Ala237_Pro161)`
fn protein_uncertain_interval(input: &str) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    delimited(char('('), protein_interval, char(')')).parse(input)
}

/// Parses protein locations written with uncertain-region syntax.
///
/// Examples: `(Ala237_Pro161)`, `(Ala237_Pro161)_(Gly170_Leu180)`
fn protein_uncertain_location(
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
fn protein_coordinate(input: &str) -> ParseResult<'_, ProteinCoordinate> {
    map(pair(protein_symbol, parse_i32), |(residue, ordinal)| {
        ProteinCoordinate { residue, ordinal }
    })
    .parse(input)
}

/// Parses protein repeat copy syntax after a protein location.
///
/// Examples: `[10]`, `[(70_80)]`
fn protein_repeat(input: &str) -> ParseResult<'_, ProteinEditKind> {
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
fn protein_extension_edit(input: &str) -> ParseResult<'_, ProteinEditKind> {
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
fn protein_n_terminal_extension_ordinal(input: &str) -> ParseResult<'_, i32> {
    map(preceded(char('-'), parse_i32), |ordinal| -ordinal).parse(input)
}

/// Parses the residue replacing the reference stop codon in C-terminal extension syntax.
///
/// Example: `Gln`
fn protein_extension_residue(input: &str) -> ParseResult<'_, String> {
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
fn protein_c_terminal_extension_state(input: &str) -> ParseResult<'_, Option<i32>> {
    preceded(
        tag("ext"),
        alt((
            value(None, pair(alt((tag("Ter"), tag("*"))), char('?'))),
            map(preceded(alt((tag("Ter"), tag("*"))), parse_i32), Some),
        )),
    )
    .parse(input)
}

/// Parses the explicitly written first residue in long protein frameshift syntax.
///
/// Example: `Pro`
fn protein_frameshift_residue(input: &str) -> ParseResult<'_, ResidueChange> {
    let (input, residue) = alt((
        map(
            delimited(char('('), protein_alternative_residues, char(')')),
            ResidueChange::Alternative,
        ),
        map(protein_symbol, ResidueChange::Known),
    ))
    .parse(input)?;

    let has_terminating_residue = match &residue {
        ResidueChange::Known(residue) => residue == "Ter",
        ResidueChange::Alternative(residues) => residues.iter().any(|residue| residue == "Ter"),
    };

    if has_terminating_residue {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )))
    } else {
        Ok((input, residue))
    }
}

/// Parses short and long protein frameshift syntax.
///
/// Examples: `fs`, `ProfsTer23`, `ProfsTer?`
fn protein_frameshift_edit(input: &str) -> ParseResult<'_, ProteinEditKind> {
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
fn protein_frameshift_stop(input: &str) -> ParseResult<'_, ProteinFrameshiftStop> {
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

fn protein_insertion_sequence(input: &str) -> ParseResult<'_, ProteinInsertionSequence> {
    alt((
        map(
            pair(
                tag("Xaa"),
                opt(delimited(
                    char('['),
                    map_res(digit1, str::parse::<usize>),
                    char(']'),
                )),
            ),
            |(_, count)| ProteinInsertionSequence::Unknown {
                count: count.unwrap_or(1),
            },
        ),
        map(
            preceded(
                alt((tag("*"), tag("Ter"))),
                map_res(digit1, str::parse::<usize>),
            ),
            |ordinal| ProteinInsertionSequence::Terminating { ordinal },
        ),
        map(protein_sequence, ProteinInsertionSequence::Known),
    ))
    .parse(input)
}

/// Parses a contiguous protein sequence.
///
/// Examples: `Ala`, `GlnSerLys`
fn protein_sequence(input: &str) -> ParseResult<'_, ProteinSequence> {
    map(many1(protein_symbol), |residues| ProteinSequence {
        residues,
    })
    .parse(input)
}

/// Parses one supported amino-acid symbol.
///
/// Examples: `Trp`, `W`, `*`
fn protein_symbol(input: &str) -> ParseResult<'_, String> {
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

fn normalize_protein_symbol(symbol: &str) -> String {
    if symbol == "*" {
        "Ter".to_string()
    } else {
        symbol.to_string()
    }
}

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;

    use super::*;

    #[test]
    fn parses_unknown_protein_insertion_sequence() {
        let (_, outcome) = all_consuming(protein_outcome)
            .parse("Arg78_Gly79insXaa[23]")
            .unwrap();

        let ProteinOutcome::Produced {
            edit: ProteinEditForm::Single(edit),
            certainty: OutcomeCertainty::Certain,
        } = outcome
        else {
            panic!("expected certain single protein edit");
        };

        assert!(matches!(
            edit.kind,
            ProteinEditKind::Insertion {
                sequence: ProteinInsertionSequence::Unknown { count: 23 },
            }
        ));
    }

    #[test]
    fn parses_bare_unknown_protein_insertion_sequence() {
        let (_, outcome) = all_consuming(protein_outcome)
            .parse("(Ser332_Ser333insXaa)")
            .unwrap();

        let ProteinOutcome::Produced {
            edit: ProteinEditForm::Single(edit),
            certainty: OutcomeCertainty::Predicted,
        } = outcome
        else {
            panic!("expected predicted single protein edit");
        };

        assert!(matches!(
            edit.kind,
            ProteinEditKind::Insertion {
                sequence: ProteinInsertionSequence::Unknown { count: 1 },
            }
        ));
    }

    #[test]
    fn parses_terminating_protein_insertion_sequence() {
        let (_, outcome) = all_consuming(protein_outcome)
            .parse("Gln746_Lys747ins*63")
            .unwrap();

        let ProteinOutcome::Produced {
            edit: ProteinEditForm::Single(edit),
            certainty: OutcomeCertainty::Certain,
        } = outcome
        else {
            panic!("expected certain single protein edit");
        };

        assert!(matches!(
            edit.kind,
            ProteinEditKind::Insertion {
                sequence: ProteinInsertionSequence::Terminating { ordinal: 63 },
            }
        ));
    }

    #[test]
    fn parses_predicted_unknown_protein_insertion() {
        let (_, outcome) = all_consuming(protein_outcome)
            .parse("(Val582_Asn583insXaa[5])")
            .unwrap();

        let ProteinOutcome::Produced {
            edit: ProteinEditForm::Single(edit),
            certainty: OutcomeCertainty::Predicted,
        } = outcome
        else {
            panic!("expected predicted single protein edit");
        };

        assert!(matches!(
            edit.kind,
            ProteinEditKind::Insertion {
                sequence: ProteinInsertionSequence::Unknown { count: 5 },
            }
        ));
    }

    #[test]
    fn parses_protein_substitution_with_alternative_residues() {
        let (_, outcome) = all_consuming(protein_outcome)
            .parse("(Gly719Ala^Ser)")
            .unwrap();

        let ProteinOutcome::Produced {
            edit: ProteinEditForm::Single(edit),
            certainty: OutcomeCertainty::Predicted,
        } = outcome
        else {
            panic!("expected predicted single protein edit");
        };

        let ProteinEditKind::Substitution {
            to: ResidueChange::Alternative(alternatives),
        } = edit.kind
        else {
            panic!("expected substitution with alternative residues");
        };

        assert_eq!(alternatives, vec!["Ala", "Ser"]);
    }

    #[test]
    fn parses_multiple_protein_substitution_alternatives() {
        let (_, outcome) = all_consuming(protein_outcome)
            .parse("(Gly56Ala^Ser^Cys)")
            .unwrap();

        let ProteinOutcome::Produced {
            edit: ProteinEditForm::Single(edit),
            certainty: OutcomeCertainty::Predicted,
        } = outcome
        else {
            panic!("expected predicted single protein edit");
        };

        let ProteinEditKind::Substitution {
            to: ResidueChange::Alternative(alternatives),
        } = edit.kind
        else {
            panic!("expected substitution with alternative residues");
        };

        assert_eq!(alternatives, vec!["Ala", "Ser", "Cys"]);
    }

    #[test]
    fn parses_frameshift_with_alternative_residues() {
        let (_, outcome) = all_consuming(protein_outcome)
            .parse("Gly719(Ala^Ser)fsTer23")
            .unwrap();

        let ProteinOutcome::Produced {
            edit: ProteinEditForm::Single(edit),
            certainty: OutcomeCertainty::Certain,
        } = outcome
        else {
            panic!("expected certain single protein edit");
        };

        let ProteinEditKind::Frameshift {
            to_residue: Some(ResidueChange::Alternative(alternatives)),
            stop,
        } = edit.kind
        else {
            panic!("expected frameshift with alternative residues");
        };

        assert_eq!(alternatives, vec!["Ala", "Ser"]);
        assert_eq!(stop.ordinal, Some(23));
    }

    #[test]
    fn rejects_frameshift_with_terminating_residue_change() {
        assert!(all_consuming(protein_outcome)
            .parse("Arg97TerfsTer23")
            .is_err());
        assert!(all_consuming(protein_outcome)
            .parse("Gly719(Ala^Ter)fsTer23")
            .is_err());
    }

    #[test]
    fn parses_alternative_protein_edits() {
        let (_, outcome) = all_consuming(protein_outcome)
            .parse("(Gly23GlufsTer7^Gly23CysfsTer26)")
            .unwrap();

        let ProteinOutcome::Produced {
            edit: ProteinEditForm::Alternative(alternatives),
            certainty: OutcomeCertainty::Predicted,
        } = outcome
        else {
            panic!("expected predicted alternative protein edits");
        };

        assert_eq!(alternatives.len(), 2);

        assert!(matches!(
            &alternatives[0].kind,
            ProteinEditKind::Frameshift {
                to_residue: Some(ResidueChange::Known(residue)),
                stop,
            } if residue == "Glu" && stop.ordinal == Some(7)
        ));

        assert!(matches!(
            &alternatives[1].kind,
            ProteinEditKind::Frameshift {
                to_residue: Some(ResidueChange::Known(residue)),
                stop,
            } if residue == "Cys" && stop.ordinal == Some(26)
        ));
    }

    #[test]
    fn parses_multiple_alternative_protein_edits() {
        let (_, form) = all_consuming(protein_edit_form)
            .parse("Gly23GlufsTer7^Gly23CysfsTer26^Gly23SerfsTer10")
            .unwrap();

        let ProteinEditForm::Alternative(alternatives) = form else {
            panic!("expected alternative protein edit form");
        };

        assert_eq!(alternatives.len(), 3);
    }

    #[test]
    fn distinguishes_residue_and_edit_alternatives() {
        let (_, residue_alternative) = all_consuming(protein_edit_form)
            .parse("Gly719Ala^Ser")
            .unwrap();

        let ProteinEditForm::Single(edit) = residue_alternative else {
            panic!("residue alternatives must remain one protein edit");
        };

        assert!(matches!(
            edit.kind,
            ProteinEditKind::Substitution {
                to: ResidueChange::Alternative(_),
            }
        ));

        let (_, edit_alternative) = all_consuming(protein_edit_form)
            .parse("Gly23GlufsTer7^Gly23CysfsTer26")
            .unwrap();

        assert!(matches!(edit_alternative, ProteinEditForm::Alternative(_)));
    }
}
