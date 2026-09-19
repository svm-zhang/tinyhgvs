//! Repeat unit, copy-number, and repeat-edit parsers.

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

/// Parses a known repeat unit.
///
/// Examples: `CAG`, `CTG`
pub(super) fn known_repeat_unit(input: &str) -> ParseResult<'_, RepeatSequenceUnit> {
    map(
        // CAG, CTG. `N` and `n` are reserved for unknown repeat units.
        verify(nucleotide_literal, |seq: &String| seq != "N" && seq != "n"),
        |seq| RepeatSequenceUnit::Known(LiteralSequenceItem { value: (seq) }),
    )
    .parse(input)
}

/// Parses an unknown repeat unit.
///
/// Examples: `N`, `n`
pub(super) fn unknown_repeat_unit(input: &str) -> ParseResult<'_, RepeatSequenceUnit> {
    // N, n
    value(RepeatSequenceUnit::Unknown, one_of("Nn")).parse(input)
}

/// Parses a known repeat copy count.
///
/// Example: `[12]`
pub(super) fn known_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
    // [12]
    delimited(
        char('['),
        map(parse_quantity, |count| Quantity::Known { count }),
        char(']'),
    )
    .parse(input)
}

/// Parses an unknown repeat copy count.
///
/// Example: `[?]`
pub(super) fn unknown_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
    // [?]
    delimited(char('['), value(Quantity::Unknown, char('?')), char(']')).parse(input)
}

/// Parses a ranged repeat copy count.
///
/// Examples: `[(60_80)]`, `[(?_60)]`, `[(60_?)]`
pub(super) fn uncertain_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
    // [(60_80)], [(?_60)], [(60_?)]
    let (input, (lo, hi)) = delimited(
        char('['),
        delimited(
            char('('),
            separated_pair(
                alt((map(parse_quantity, Some), value(None, char('?')))),
                char('_'),
                alt((map(parse_quantity, Some), value(None, char('?')))),
            ),
            char(')'),
        ),
        char(']'),
    )
    .parse(input)?;

    if lo.is_none() && hi.is_none() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    Ok((input, Quantity::Uncertain { lo, hi }))
}

/// Parses a repeat edit with a known repeat unit.
///
/// Examples: `CAG[23]`, `CAG[(60_80)]`, `CAG[?]`
pub(super) fn known_repeat_edit(input: &str) -> ParseResult<'_, RepeatEdit> {
    // CAG[23], CAG[(60_80)], CAG[?]
    map(
        pair(
            known_repeat_unit,
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

/// Parses a repeat edit with an unknown repeat unit.
///
/// Examples: `N[12]`, `N[(60_80)]`, `N[?]`
pub(super) fn unknown_repeat_edit(input: &str) -> ParseResult<'_, RepeatEdit> {
    // N[12], N[(60_80)], N[?]
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

/// Parses top-level nucleotide repeat edits.
///
/// Examples: `CTG[9]TTG[1]CTG[13]`, `[14]`, `[(80_100)]`
pub(super) fn repeat_edits(input: &str) -> ParseResult<'_, NucleotideEditKind> {
    map(
        alt((
            // CTG[9]TTG[1]CTG[13]
            many1(known_repeat_edit),
            // [14], [(80_100)]
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
