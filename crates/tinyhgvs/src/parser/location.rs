//! Nucleotide location and coordinate parsers.

use nom::branch::alt;
use nom::character::complete::char;
use nom::combinator::verify;
use nom::combinator::{map, opt, value};
use nom::sequence::{delimited, pair, preceded};
use nom::Parser;

use super::core::{parse_i32, parse_position, range_with};
use super::ParseResult;
use crate::model::{Interval, Location, NucleotideAnchor, NucleotideCoordinate};

/// Parses a nucleotide location as known or uncertain.
///
/// Examples: `93`, `93_94`, `(71_72)`, `(123_234)_(345_456)`
pub(super) fn nucleotide_location(input: &str) -> ParseResult<'_, Location<NucleotideCoordinate>> {
    // Reject location description such as `(?_?)`, `(?_?)_(?_?)`.
    let is_valid_uncertain_location = |loc: &Interval<Interval<NucleotideCoordinate>>| {
        !(loc.start.is_fully_unknown() && loc.end.as_ref().map_or(true, Interval::is_fully_unknown))
    };
    // Reject partially unknown known-location intervals such as `?_87` and
    // `123_?`, while keeping whole-location `?_?`.
    let is_valid_known_interval = |interval: &Interval<NucleotideCoordinate>| {
        !interval.has_unknown_bound() || interval.is_fully_unknown()
    };

    alt((
        // (71_72) and (123_234)_(345_456), (?_87), (123_?)_(?_456)
        map(
            verify(nucleotide_uncertain_location, is_valid_uncertain_location),
            Location::from_uncertain,
        ),
        // 93 and 93_94, plus whole-location `?_?`
        map(
            verify(nucleotide_interval, is_valid_known_interval),
            Location::from_known,
        ),
    ))
    .parse(input)
}

/// Parses nucleotide location as a single position/coordinate or an interval
/// joined by `_`.
/// - `93`
/// - `93_94`
/// - `?_87`
/// - `123_?`
/// - `?_?`
pub(super) fn nucleotide_interval(input: &str) -> ParseResult<'_, Interval<NucleotideCoordinate>> {
    alt((
        // Interval coordinate, `93_94`
        |input| range_with(input, nucleotide_coordinate),
        // Single position coordinate, `93`
        map(nucleotide_coordinate, |start| Interval { start, end: None }),
    ))
    .parse(input)
}

/// Parses one uncertain interval unit with parentheses.
///
/// - `(71_72)`
/// - `(123_?)`
/// - `(?_87)`
pub(super) fn nucleotide_uncertain_interval(
    input: &str,
) -> ParseResult<'_, Interval<NucleotideCoordinate>> {
    delimited(char('('), nucleotide_interval, char(')')).parse(input)
}

/// Parses a nucleotide uncertain location.
///
/// One location can be either one or two uncertain interval units separated by
/// `_`.
pub(super) fn nucleotide_uncertain_location(
    input: &str,
) -> ParseResult<'_, Interval<Interval<NucleotideCoordinate>>> {
    alt((
        // (123_234)_(345_456)
        |input| range_with(input, nucleotide_uncertain_interval),
        // (71_72)
        map(nucleotide_uncertain_interval, |start| Interval {
            start,
            end: None,
        }),
    ))
    .parse(input)
}

/// Parses a nucleotide coordinate with anchor and optional offset.
fn nucleotide_coordinate(input: &str) -> ParseResult<'_, NucleotideCoordinate> {
    alt((
        // -18, -106+2, -84-1
        map(
            pair(
                preceded(char('-'), parse_position),
                opt(pair(alt((char('+'), char('-'))), parse_i32)),
            ),
            |(coordinate, offset)| {
                let offset = offset
                    .map(|(sign, value)| if sign == '-' { -value } else { value })
                    .unwrap_or(0);

                NucleotideCoordinate::known(NucleotideAnchor::RelativeCdsStart, -coordinate, offset)
            },
        ),
        // *18, *639-1
        map(
            pair(
                preceded(char('*'), parse_position),
                opt(pair(alt((char('+'), char('-'))), parse_i32)),
            ),
            |(coordinate, offset)| {
                let offset = offset
                    .map(|(sign, value)| if sign == '-' { -value } else { value })
                    .unwrap_or(0);

                NucleotideCoordinate::known(NucleotideAnchor::RelativeCdsEnd, coordinate, offset)
            },
        ),
        // 93, 93+1, 93-2
        map(
            pair(parse_i32, opt(pair(alt((char('+'), char('-'))), parse_i32))),
            |(coordinate, offset)| {
                let offset = offset
                    .map(|(sign, value)| if sign == '-' { -value } else { value })
                    .unwrap_or(0);

                NucleotideCoordinate::known(NucleotideAnchor::Absolute, coordinate, offset)
            },
        ),
        // ?, unknown coordinate
        value(NucleotideCoordinate::Unknown, char('?')),
    ))
    .parse(input)
}

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;
    use nom::Parser;

    use super::*;

    #[test]
    fn parses_nucleotide_position_branches() {
        let (_, coding) = all_consuming(nucleotide_coordinate).parse("93+1").unwrap();
        assert_eq!(coding.anchor().unwrap(), NucleotideAnchor::Absolute);
        assert_eq!(coding.coordinate(), Some(93));
        assert_eq!(coding.offset().unwrap(), 1);

        let (_, upstream_intronic) = all_consuming(nucleotide_coordinate).parse("93-2").unwrap();
        assert_eq!(
            upstream_intronic.anchor().unwrap(),
            NucleotideAnchor::Absolute
        );
        assert_eq!(upstream_intronic.coordinate(), Some(93));
        assert_eq!(upstream_intronic.offset().unwrap(), -2);

        let (_, utr5) = all_consuming(nucleotide_coordinate).parse("-18").unwrap();
        assert_eq!(utr5.anchor().unwrap(), NucleotideAnchor::RelativeCdsStart);
        assert_eq!(utr5.coordinate(), Some(-18));
        assert_eq!(utr5.offset().unwrap(), 0);

        let (_, utr5_intronic) = all_consuming(nucleotide_coordinate)
            .parse("-106+2")
            .unwrap();
        assert_eq!(
            utr5_intronic.anchor().unwrap(),
            NucleotideAnchor::RelativeCdsStart
        );
        assert_eq!(utr5_intronic.coordinate(), Some(-106));
        assert_eq!(utr5_intronic.offset().unwrap(), 2);

        let (_, utr5_intronic_upstream) =
            all_consuming(nucleotide_coordinate).parse("-84-1").unwrap();
        assert_eq!(
            utr5_intronic_upstream.anchor().unwrap(),
            NucleotideAnchor::RelativeCdsStart
        );
        assert_eq!(utr5_intronic_upstream.coordinate().unwrap(), -84);
        assert_eq!(utr5_intronic_upstream.offset().unwrap(), -1);

        let (_, utr3) = all_consuming(nucleotide_coordinate).parse("*18").unwrap();
        assert_eq!(utr3.anchor().unwrap(), NucleotideAnchor::RelativeCdsEnd);
        assert_eq!(utr3.coordinate(), Some(18));
        assert_eq!(utr3.offset().unwrap(), 0);

        let (_, utr3_intronic) = all_consuming(nucleotide_coordinate)
            .parse("*639-1")
            .unwrap();
        assert_eq!(
            utr3_intronic.anchor().unwrap(),
            NucleotideAnchor::RelativeCdsEnd
        );
        assert_eq!(utr3_intronic.coordinate(), Some(639));
        assert_eq!(utr3_intronic.offset().unwrap(), -1);

        let (_, unknown) = all_consuming(nucleotide_coordinate).parse("?").unwrap();
        assert!(unknown.is_unknown());
        assert_eq!(unknown.anchor(), None);
        assert_eq!(unknown.coordinate(), None);
        assert_eq!(unknown.offset(), None);

        assert!(all_consuming(nucleotide_coordinate).parse("-0").is_err());
        assert!(all_consuming(nucleotide_coordinate).parse("-0+2").is_err());
        assert!(all_consuming(nucleotide_coordinate).parse("*0").is_err());
        assert!(all_consuming(nucleotide_coordinate).parse("*0-1").is_err());
    }
}
