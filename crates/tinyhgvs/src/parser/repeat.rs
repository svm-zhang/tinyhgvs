use nom::branch::alt;
use nom::character::complete::{char, one_of};
use nom::combinator::{map, value, verify};
use nom::multi::many1;
use nom::sequence::{delimited, pair, separated_pair};
use nom::Parser;

use super::core::{nucleotide_literal, parse_quantity};
use super::ParseResult;

use crate::model::{
    LiteralSequenceItem, NucleotideEditKind, Quantity, RepeatEdit, RepeatSequenceUnit,
};

pub(super) fn known_repeat_unit(input: &str) -> ParseResult<'_, RepeatSequenceUnit> {
    map(
        // Because nucleotide_literal consumes any alphabetic string, it also
        // consumes "N" or "n", the verify function makes sure when that happens,
        // this parser fails. This makes sure `N[12]` is captured by the
        // unknown_repeat_edit correctly. Without using the verify function,
        //  I need to swap the order inside alt as: alt((unknown, known))
        //  inside nucleotide_sequence_item parser.
        verify(nucleotide_literal, |seq: &String| seq != "N" && seq != "n"),
        |seq| RepeatSequenceUnit::Known(LiteralSequenceItem { value: (seq) }),
    )
    .parse(input)
}

pub(super) fn unknown_repeat_unit(input: &str) -> ParseResult<'_, RepeatSequenceUnit> {
    value(RepeatSequenceUnit::Unknown, one_of("Nn")).parse(input)
}

pub(super) fn known_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
    delimited(
        char('['),
        map(parse_quantity, |count| Quantity::Known { count }),
        char(']'),
    )
    .parse(input)
}

pub(super) fn unknown_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
    delimited(char('['), value(Quantity::Unknown, char('?')), char(']')).parse(input)
}

pub(super) fn uncertain_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
    delimited(
        char('['),
        // (A_B), (A_?), (?_B), (?_?)
        delimited(
            char('('),
            map(
                separated_pair(
                    alt((map(parse_quantity, Some), value(None, char('?')))),
                    char('_'),
                    alt((map(parse_quantity, Some), value(None, char('?')))),
                ),
                |(lo, hi)| Quantity::Uncertain { lo, hi },
            ),
            char(')'),
        ),
        char(']'),
    )
    .parse(input)
}

pub(super) fn known_repeat_edit(input: &str) -> ParseResult<'_, RepeatEdit> {
    map(
        pair(
            known_repeat_unit,
            alt((known_repeat_copy, uncertain_repeat_copy)),
        ),
        |(unit, quantity)| RepeatEdit {
            quantity,
            unit: Some(unit),
        },
    )
    .parse(input)
}

pub(super) fn unknown_repeat_edit(input: &str) -> ParseResult<'_, RepeatEdit> {
    map(
        pair(
            unknown_repeat_unit,
            alt((
                known_repeat_copy,
                uncertain_repeat_copy,
                unknown_repeat_copy,
            )),
        ),
        |(unit, quantity)| RepeatEdit {
            quantity,
            unit: Some(unit),
        },
    )
    .parse(input)
}

// Note: I do not need to add the unknown_repeat_edit pattern as those
// are mainly present as the inserted or replaced item in ins and delins.
// unknown_repeat_edit parses N[80], N[(80_100)], and N[?]
pub(super) fn repeat_edits(input: &str) -> ParseResult<'_, NucleotideEditKind> {
    map(
        alt((
            // CTG[9]TTG[1]CTG[13]
            many1(known_repeat_edit),
            // one occurrence of either [14] or [(80_100)]
            map(
                alt((uncertain_repeat_copy, known_repeat_copy)),
                |quantity| {
                    vec![RepeatEdit {
                        unit: None,
                        quantity,
                    }]
                },
            ),
        )),
        |blocks| NucleotideEditKind::Repeat { blocks },
    )
    .parse(input)
}
