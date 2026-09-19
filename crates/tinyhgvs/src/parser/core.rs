//! Shared grammar primitives used by multiple parser families.

use nom::branch::alt;
use nom::bytes::complete::take_while1;
use nom::character::complete::{char, digit1};
use nom::combinator::{map, map_res, opt, value};
use nom::sequence::{delimited, pair, separated_pair};
use nom::Parser;

use super::ParseResult;
use crate::model::{Accession, CoordinateSystem, Interval, ReferenceSpec};

/// Parses the HGVS reference identifier field into a model::ReferenceSpec type.
/// Genomic reference plus a transcript context form is supported.
pub(super) fn reference_spec(input: &str) -> ParseResult<'_, ReferenceSpec> {
    map(
        pair(accession, opt(delimited(char('('), accession, char(')')))),
        |(primary, context)| ReferenceSpec {
            primary: Accession::new(primary),
            context: context.map(Accession::new),
        },
    )
    .parse(input)
}

/// Parses sequence accession such as `NM_004006.2` or `ENST00000351052.5`.
fn accession(input: &str) -> ParseResult<'_, String> {
    map(
        take_while1(|c: char| c.is_ascii_alphanumeric() || matches!(c, '_' | '.')),
        str::to_string,
    )
    .parse(input)
}

/// Parses the one-letter HGVS coordinate system marker.
pub(super) fn coordinate_system(input: &str) -> ParseResult<'_, CoordinateSystem> {
    alt((
        value(CoordinateSystem::Genomic, char('g')),
        value(CoordinateSystem::CodingDna, char('c')),
        value(CoordinateSystem::Rna, char('r')),
        value(CoordinateSystem::Protein, char('p')),
    ))
    .parse(input)
}

/// Parses a reusable range surface written as `thing_thing`.
pub(super) fn range_with<T, P>(input: &str, parse_item: P) -> ParseResult<'_, Interval<T>>
where
    P: Copy + Fn(&str) -> ParseResult<'_, T>,
{
    map(
        separated_pair(parse_item, char('_'), parse_item),
        |(start, end)| Interval {
            start,
            end: Some(end),
        },
    )
    .parse(input)
}

/// Parses an unsigned decimal integer into `i32`.
pub(super) fn parse_i32(input: &str) -> ParseResult<'_, i32> {
    map_res(digit1, str::parse::<i32>).parse(input)
}

/// Parses a one-based coordinate position.
///
/// This rejects `0`, which is not a valid HGVS position.
pub(super) fn parse_position(input: &str) -> ParseResult<'_, i32> {
    let (input, value) = parse_i32(input)?;
    if value == 0 {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )))
    } else {
        Ok((input, value))
    }
}

/// Parses an unsigned decimal integer for counts and amounts.
pub(super) fn parse_quantity(input: &str) -> ParseResult<'_, usize> {
    map_res(digit1, str::parse::<usize>).parse(input)
}

/// Parses a literal nucleotide token.
///
/// Examples: `G`, `AGGG`, `CAG`
pub(super) fn nucleotide_literal(input: &str) -> ParseResult<'_, String> {
    map(
        take_while1(|c: char| c.is_ascii_alphabetic()),
        str::to_string,
    )
    .parse(input)
}
