//! HGVS variant parsers.
//!
//! The parsing strategy:
//!
//! - parse the accepted syntax into the Rust data model
//! - route rejected inputs to the lightweight diagnostic classifier

use nom::branch::alt;
use nom::bytes::complete::{tag, take_while1};
use nom::character::complete::{char, digit1};
use nom::combinator::{all_consuming, map, map_res, opt, value};
use nom::multi::{many1, separated_list1};
use nom::sequence::{delimited, pair, preceded};
use nom::{IResult, Parser};

use crate::diagnostics::classify_parse_failure;
use crate::error::ParseHgvsError;
use crate::model::{
    Accession, CoordinateSystem, CopiedSequenceItem, HgvsVariant, Interval, LiteralSequenceItem,
    NucleotideAnchor, NucleotideCoordinate, NucleotideEdit, NucleotideSequenceItem,
    NucleotideVariant, ProteinCoordinate, ProteinEdit, ProteinEffect, ProteinSequence,
    ProteinVariant, ReferenceSpec, RepeatSequenceItem, VariantDescription,
};

type ParseResult<'a, T> = IResult<&'a str, T>;

const PROTEIN_SYMBOLS: &[&str] = &[
    "Ter", "Sec", "Pyl", "Xaa", "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His",
    "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "*", "A", "R",
    "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
];

/// Parses an HGVS string into the Rust [`HgvsVariant`] model.
///
/// Leading and trailing whitespace are ignored before parsing.
///
/// The returned model keeps the HGVS expression split into:
///
/// - `reference`: the reference source for a variant.
/// - `coordinate_system`: the one-letter reference coordinate type.
/// - `description`, the nucleotide or protein variant description.
///
/// # Examples
///
/// A splice-adjacent substitution in an intron:
///
/// ```rust
/// use tinyhgvs::{NucleotideAnchor, NucleotideEdit, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("  NM_004006.2:c.357+1G>A  ").unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         assert_eq!(nucleotide.location.start.anchor, NucleotideAnchor::Absolute);
///         assert_eq!(nucleotide.location.start.coordinate, 357);
///         assert_eq!(nucleotide.location.start.offset, 1);
///         assert!(matches!(
///             nucleotide.edit,
///             NucleotideEdit::Substitution { ref reference, ref alternate }
///                 if reference == "G" && alternate == "A"
///         ));
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
///
/// A 5' UTR substitution keeps the signed coordinate from the HGVS string:
///
/// ```rust
/// use tinyhgvs::{NucleotideAnchor, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_007373.4:c.-1C>T").unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         assert_eq!(nucleotide.location.start.anchor, NucleotideAnchor::RelativeCdsStart);
///         assert_eq!(nucleotide.location.start.coordinate, -1);
///         assert_eq!(nucleotide.location.start.offset, 0);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
///
/// A nonsense mutation leading to an early termination consequence at protein-level:
///
/// ```rust
/// use tinyhgvs::{ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.1:p.Trp24Ter").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(protein) => {
///         assert!(!protein.is_predicted);
///         assert!(matches!(protein.effect, ProteinEffect::Edit { .. }));
///     }
///     _ => unreachable!("expected protein variant"),
/// }
/// ```
///
/// Unsupported syntax is reported as a structured [`crate::ParseHgvsError`]:
///
/// ```rust
/// use tinyhgvs::parse_hgvs;
///
/// let error = parse_hgvs("NM_004006.3:r.spl").unwrap_err();
/// assert_eq!(error.code(), "unsupported.rna_special_state");
/// ```
pub fn parse_hgvs(input: &str) -> Result<HgvsVariant, ParseHgvsError> {
    // Trim leading and trailing spaces.
    let input = input.trim();
    all_consuming(hgvs_variant)
        .parse(input)
        .map(|(_, variant)| variant)
        .map_err(|_| classify_parse_failure(input))
}

/// Parses either a variant with reference identifier or not. Context-dependent
/// shorthand protein-level description is allowed, e.g. "p.Gly12Asp".
fn hgvs_variant(input: &str) -> ParseResult<'_, HgvsVariant> {
    // Match either a full nucleotide or shorthand protein syntax.
    alt((variant_with_reference, protein_variant_without_reference)).parse(input)
}

/// Parses the full HGVS variant (with a reference identifier).
fn variant_with_reference(input: &str) -> ParseResult<'_, HgvsVariant> {
    // Parses the reference field
    let (input, reference) = reference_spec(input)?;
    // Reads the separator between reference and coordinate type
    let (input, _) = char(':')(input)?;
    // Parses the coordinate type, e.g. "g", "c", "r", etc
    let (input, coordinate_system) = coordinate_system(input)?;
    // Reads the separator to move into description
    let (input, _) = char('.')(input)?;
    // Parses the description syntax
    let (input, description) = variant_description(coordinate_system, input)?;

    Ok((
        input,
        build_variant(Some(reference), coordinate_system, description),
    ))
}

/// Parses context-dependent shorthand protein-level variant.
fn protein_variant_without_reference(input: &str) -> ParseResult<'_, HgvsVariant> {
    let (input, _) = tag("p.")(input)?;
    let (input, description) = protein_description(input)?;
    Ok((
        input,
        build_variant(None, CoordinateSystem::Protein, description),
    ))
}

/// Builds an [`HgvsVariant`] from its structural components.
fn build_variant(
    reference: Option<ReferenceSpec>,
    coordinate_system: CoordinateSystem,
    description: VariantDescription,
) -> HgvsVariant {
    HgvsVariant {
        reference,
        coordinate_system,
        description,
    }
}

/// Route a variant to its matching parser based on coordinate type.
fn variant_description(
    coordinate_system: CoordinateSystem,
    input: &str,
) -> ParseResult<'_, VariantDescription> {
    if coordinate_system.is_protein() {
        protein_description(input)
    } else {
        nucleotide_description(input)
    }
}

/// Parses the HGVS reference identifier field into a model::ReferenceSpec type.
/// Genomic reference plus a transcript context form is supported.
fn reference_spec(input: &str) -> ParseResult<'_, ReferenceSpec> {
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
fn coordinate_system(input: &str) -> ParseResult<'_, CoordinateSystem> {
    alt((
        value(CoordinateSystem::Genomic, char('g')),
        value(CoordinateSystem::CircularGenomic, char('o')),
        value(CoordinateSystem::Mitochondrial, char('m')),
        value(CoordinateSystem::CodingDna, char('c')),
        value(CoordinateSystem::NonCodingDna, char('n')),
        value(CoordinateSystem::Rna, char('r')),
        value(CoordinateSystem::Protein, char('p')),
    ))
    .parse(input)
}

/// Parses nucleotide description composed of location plus edit.
fn nucleotide_description(input: &str) -> ParseResult<'_, VariantDescription> {
    map(
        pair(nucleotide_interval, nucleotide_edit),
        |(location, edit)| VariantDescription::Nucleotide(NucleotideVariant { location, edit }),
    )
    .parse(input)
}

/// Parses nucleotide location as a single position/coordinate or an interval
/// joined by `_`.
fn nucleotide_interval(input: &str) -> ParseResult<'_, Interval<NucleotideCoordinate>> {
    alt((
        map(
            pair(
                nucleotide_coordinate,
                preceded(char('_'), nucleotide_coordinate),
            ),
            |(start, end)| Interval {
                start,
                end: Some(end),
            },
        ),
        map(nucleotide_coordinate, |start| Interval { start, end: None }),
    ))
    .parse(input)
}

/// Parses a nucleotide coordinate with anchor and optional offset.
fn nucleotide_coordinate(input: &str) -> ParseResult<'_, NucleotideCoordinate> {
    alt((
        map(preceded(char('-'), parse_i32), |coordinate| {
            NucleotideCoordinate {
                anchor: NucleotideAnchor::RelativeCdsStart,
                coordinate: -coordinate,
                offset: 0,
            }
        }),
        map(preceded(char('*'), parse_i32), |coordinate| {
            NucleotideCoordinate {
                anchor: NucleotideAnchor::RelativeCdsEnd,
                coordinate,
                offset: 0,
            }
        }),
        map(
            pair(parse_i32, opt(pair(alt((char('+'), char('-'))), parse_i32))),
            |(coordinate, offset)| {
                let offset = offset
                    .map(|(sign, value)| if sign == '-' { -value } else { value })
                    .unwrap_or(0);

                NucleotideCoordinate {
                    anchor: NucleotideAnchor::Absolute,
                    coordinate,
                    offset,
                }
            },
        ),
    ))
    .parse(input)
}

/// Parses an unsigned decimal integer into `i32`.
fn parse_i32(input: &str) -> ParseResult<'_, i32> {
    map_res(digit1, str::parse::<i32>).parse(input)
}

/// Parses an unsigned decimal integer into `usize`.
fn parse_usize(input: &str) -> ParseResult<'_, usize> {
    map_res(digit1, str::parse::<usize>).parse(input)
}

/// Parses the currently supported nucleotide edit families.
fn nucleotide_edit(input: &str) -> ParseResult<'_, NucleotideEdit> {
    alt((
        value(NucleotideEdit::NoChange, char('=')),
        map(
            preceded(tag("delins"), nucleotide_sequence_items),
            |items| NucleotideEdit::DeletionInsertion { items },
        ),
        value(NucleotideEdit::Deletion, tag("del")),
        value(NucleotideEdit::Duplication, tag("dup")),
        map(preceded(tag("ins"), nucleotide_sequence_items), |items| {
            NucleotideEdit::Insertion { items }
        }),
        value(NucleotideEdit::Inversion, tag("inv")),
        map(
            pair(nucleotide_literal, preceded(char('>'), nucleotide_literal)),
            |(reference, alternate)| NucleotideEdit::Substitution {
                reference,
                alternate,
            },
        ),
    ))
    .parse(input)
}

/// Parses inserted or replacement sequence items in an `ins` or `delins` variant.
fn nucleotide_sequence_items(input: &str) -> ParseResult<'_, Vec<NucleotideSequenceItem>> {
    map(
        alt((
            delimited(
                char('['),
                separated_list1(char(';'), nucleotide_sequence_item),
                char(']'),
            ),
            map(nucleotide_sequence_item, |item| vec![item]),
        )),
        |items| items,
    )
    .parse(input)
}

/// Parses one sequence item as literal, repeat, or copied sequence.
fn nucleotide_sequence_item(input: &str) -> ParseResult<'_, NucleotideSequenceItem> {
    alt((
        map(sequence_repeat, NucleotideSequenceItem::Repeat),
        map(sequence_segment, NucleotideSequenceItem::Copied),
        map(nucleotide_literal, |value| {
            NucleotideSequenceItem::Literal(LiteralSequenceItem { value })
        }),
    ))
    .parse(input)
}

/// Parses repeat expression such as `T[12]`.
fn sequence_repeat(input: &str) -> ParseResult<'_, RepeatSequenceItem> {
    map(
        pair(
            nucleotide_literal,
            delimited(char('['), parse_usize, char(']')),
        ),
        |(unit, count)| RepeatSequenceItem { unit, count },
    )
    .parse(input)
}

/// Parses a segment- or interval-type edit component that comes from either
/// local (current) or remote (other) reference source.
fn sequence_segment(input: &str) -> ParseResult<'_, CopiedSequenceItem> {
    alt((remote_sequence_segment, current_reference_sequence_segment)).parse(input)
}

/// Parses a current-reference segment such as `850_900inv`.
fn current_reference_sequence_segment(input: &str) -> ParseResult<'_, CopiedSequenceItem> {
    map(
        pair(nucleotide_interval, opt(tag("inv"))),
        |(source_location, is_inverted)| CopiedSequenceItem {
            source_reference: None,
            source_coordinate_system: None,
            source_location,
            is_inverted: is_inverted.is_some(),
        },
    )
    .parse(input)
}

/// Parses a other-reference segment such as `NC_000022.10:g.35788169_35788352`.
fn remote_sequence_segment(input: &str) -> ParseResult<'_, CopiedSequenceItem> {
    map(
        (
            reference_spec,
            char(':'),
            coordinate_system,
            char('.'),
            nucleotide_interval,
            opt(tag("inv")),
        ),
        |(source_reference, _, source_coordinate_system, _, source_location, is_inverted)| {
            CopiedSequenceItem {
                source_reference: Some(source_reference),
                source_coordinate_system: Some(source_coordinate_system),
                source_location,
                is_inverted: is_inverted.is_some(),
            }
        },
    )
    .parse(input)
}

/// Parses an edit component of literal nucleotide base changes.
fn nucleotide_literal(input: &str) -> ParseResult<'_, String> {
    map(
        take_while1(|c: char| c.is_ascii_alphabetic()),
        str::to_string,
    )
    .parse(input)
}

/// Parses protein description of either known or predicted effect and consequences.
fn protein_description(input: &str) -> ParseResult<'_, VariantDescription> {
    alt((
        map(delimited(char('('), protein_effect, char(')')), |effect| {
            VariantDescription::Protein(ProteinVariant {
                is_predicted: true,
                effect,
            })
        }),
        map(protein_effect, |effect| {
            VariantDescription::Protein(ProteinVariant {
                is_predicted: false,
                effect,
            })
        }),
    ))
    .parse(input)
}

/// Parses the supported protein effect types. Differentiating unknown consequence
/// from the case where mutation produces no protein.
fn protein_effect(input: &str) -> ParseResult<'_, ProteinEffect> {
    alt((
        value(ProteinEffect::Unknown, char('?')),
        value(ProteinEffect::NoProteinProduced, char('0')),
        map(pair(protein_interval, protein_edit), |(location, edit)| {
            ProteinEffect::Edit { location, edit }
        }),
    ))
    .parse(input)
}

/// Parses a single protein position or an interval.
fn protein_interval(input: &str) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    alt((
        map(
            pair(protein_coordinate, preceded(char('_'), protein_coordinate)),
            |(start, end)| Interval {
                start,
                end: Some(end),
            },
        ),
        map(protein_coordinate, |start| Interval { start, end: None }),
    ))
    .parse(input)
}

/// Parses a protein symbol followed by its ordinal.
fn protein_coordinate(input: &str) -> ParseResult<'_, ProteinCoordinate> {
    map(pair(protein_symbol, parse_i32), |(residue, ordinal)| {
        ProteinCoordinate { residue, ordinal }
    })
    .parse(input)
}

/// Parses the currently supported protein edit families.
fn protein_edit(input: &str) -> ParseResult<'_, ProteinEdit> {
    alt((
        value(ProteinEdit::Unknown, char('?')),
        value(ProteinEdit::NoChange, char('=')),
        map(preceded(tag("delins"), protein_sequence), |sequence| {
            ProteinEdit::DeletionInsertion { sequence }
        }),
        value(ProteinEdit::Deletion, tag("del")),
        value(ProteinEdit::Duplication, tag("dup")),
        map(preceded(tag("ins"), protein_sequence), |sequence| {
            ProteinEdit::Insertion { sequence }
        }),
        map(protein_symbol, |to| ProteinEdit::Substitution { to }),
    ))
    .parse(input)
}

/// Parses a contiguous protein sequence.
fn protein_sequence(input: &str) -> ParseResult<'_, ProteinSequence> {
    map(many1(protein_symbol), |residues| ProteinSequence {
        residues,
    })
    .parse(input)
}

/// Parses one supported amino-acid symbol.
fn protein_symbol(input: &str) -> ParseResult<'_, String> {
    for symbol in PROTEIN_SYMBOLS {
        if let Some(rest) = input.strip_prefix(symbol) {
            return Ok((rest, (*symbol).to_string()));
        }
    }

    Err(nom::Err::Error(nom::error::Error::new(
        input,
        nom::error::ErrorKind::Tag,
    )))
}

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;

    use super::*;

    #[test]
    fn parses_nucleotide_position_branches() {
        let (_, coding) = all_consuming(nucleotide_coordinate).parse("93+1").unwrap();
        assert_eq!(coding.anchor, NucleotideAnchor::Absolute);
        assert_eq!(coding.coordinate, 93);
        assert_eq!(coding.offset, 1);

        let (_, upstream_intronic) = all_consuming(nucleotide_coordinate).parse("93-2").unwrap();
        assert_eq!(upstream_intronic.anchor, NucleotideAnchor::Absolute);
        assert_eq!(upstream_intronic.coordinate, 93);
        assert_eq!(upstream_intronic.offset, -2);

        let (_, utr5) = all_consuming(nucleotide_coordinate).parse("-18").unwrap();
        assert_eq!(utr5.anchor, NucleotideAnchor::RelativeCdsStart);
        assert_eq!(utr5.coordinate, -18);
        assert_eq!(utr5.offset, 0);

        let (_, utr3) = all_consuming(nucleotide_coordinate).parse("*18").unwrap();
        assert_eq!(utr3.anchor, NucleotideAnchor::RelativeCdsEnd);
        assert_eq!(utr3.coordinate, 18);
        assert_eq!(utr3.offset, 0);
    }

    #[test]
    fn parses_nucleotide_edit_branches() {
        assert_eq!(
            all_consuming(nucleotide_edit).parse("=").unwrap().1,
            NucleotideEdit::NoChange
        );
        assert_eq!(
            all_consuming(nucleotide_edit).parse("del").unwrap().1,
            NucleotideEdit::Deletion
        );
        assert!(all_consuming(nucleotide_edit).parse("delA").is_err());
        assert_eq!(
            all_consuming(nucleotide_edit).parse("dup").unwrap().1,
            NucleotideEdit::Duplication
        );
        assert_eq!(
            all_consuming(nucleotide_edit).parse("inv").unwrap().1,
            NucleotideEdit::Inversion
        );
        assert!(matches!(
            all_consuming(nucleotide_edit).parse("C>A").unwrap().1,
            NucleotideEdit::Substitution { .. }
        ));
        assert!(matches!(
            all_consuming(nucleotide_edit).parse("insT").unwrap().1,
            NucleotideEdit::Insertion { .. }
        ));
        assert!(matches!(
            all_consuming(nucleotide_edit).parse("delinsT").unwrap().1,
            NucleotideEdit::DeletionInsertion { .. }
        ));
    }

    #[test]
    fn parses_nucleotide_sequence_items() {
        let (_, literal) = all_consuming(nucleotide_sequence_items).parse("T").unwrap();
        assert_eq!(literal.len(), 1);

        let (_, repeat) = all_consuming(nucleotide_sequence_items)
            .parse("N[12]")
            .unwrap();
        assert!(matches!(
            repeat.first().unwrap(),
            NucleotideSequenceItem::Repeat(RepeatSequenceItem { unit, count })
                if unit == "N" && *count == 12
        ));

        let (_, local) = all_consuming(nucleotide_sequence_items)
            .parse("850_900inv")
            .unwrap();
        assert!(matches!(
            local.first().unwrap(),
            NucleotideSequenceItem::Copied(CopiedSequenceItem {
                source_reference: None,
                source_coordinate_system: None,
                is_inverted: true,
                ..
            })
        ));

        let (_, remote) = all_consuming(nucleotide_sequence_items)
            .parse("[NC_000022.10:g.35788169_35788352]")
            .unwrap();
        assert!(matches!(
            remote.first().unwrap(),
            NucleotideSequenceItem::Copied(CopiedSequenceItem {
                source_reference: Some(_),
                source_coordinate_system: Some(CoordinateSystem::Genomic),
                ..
            })
        ));
    }

    #[test]
    fn parses_protein_effect_branches() {
        assert_eq!(
            all_consuming(protein_effect).parse("?").unwrap().1,
            ProteinEffect::Unknown
        );
        assert_eq!(
            all_consuming(protein_effect).parse("0").unwrap().1,
            ProteinEffect::NoProteinProduced
        );
        assert!(matches!(
            all_consuming(protein_effect).parse("Met1?").unwrap().1,
            ProteinEffect::Edit {
                edit: ProteinEdit::Unknown,
                ..
            }
        ));
        assert!(matches!(
            all_consuming(protein_effect).parse("Trp24Ter").unwrap().1,
            ProteinEffect::Edit {
                edit: ProteinEdit::Substitution { .. },
                ..
            }
        ));
    }
}
