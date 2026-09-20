//! Genomic description, nucleotide edit, and nucleotide allele parsers.

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::char;
use nom::combinator::{map, opt, value};
use nom::multi::separated_list1;
use nom::sequence::{delimited, pair, preceded, separated_pair};
use nom::Parser;

use super::ParseResult;
use crate::model::{
    Allele, AlleleForm, AllelePhase, AlleleVariant, CopiedSequenceItem, GenomicOutcome,
    LiteralSequenceItem, NucleotideEdit, NucleotideEditKind, NucleotideSequenceItem,
    VariantDescription,
};

use super::core::{coordinate_system, nucleotide_literal, reference_spec};
use super::location::{nucleotide_interval, nucleotide_location};
use super::repeat::*;

/// Parses one or more nucleotide edits written inside one allele.
///
/// Example: `[123G>A;345del]`
pub(super) fn nucleotide_variants_on_allele(input: &str) -> ParseResult<'_, Vec<NucleotideEdit>> {
    delimited(
        char('['),
        separated_list1(char(';'), nucleotide_edit),
        char(']'),
    )
    .parse(input)
}

/// Parses nucleotide cis allele.
///
/// Example: `[123G>A;345del]`
fn nucleotide_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<NucleotideEdit>> {
    map(nucleotide_variants_on_allele, |variants| AlleleVariant {
        allele_one: Allele::from_variants(variants),
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

/// Parses two nucleotide alleles in trans, with optional unphased edits.
///
/// Examples: `[123G>A];[345del]`, `[A];[B](;)C`
fn nucleotide_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<NucleotideEdit>> {
    map(
        pair(
            separated_pair(
                nucleotide_variants_on_allele,
                char(';'),
                nucleotide_variants_on_allele,
            ),
            opt(preceded(
                tag("(;)"),
                separated_list1(tag("(;)"), nucleotide_edit),
            )),
        ),
        |((a1, a2), unphased)| AlleleVariant {
            allele_one: Allele::from_variants(a1),
            allele_two: Some(Allele::from_variants(a2)),
            phase: Some(AllelePhase::Trans),
            variants_unphased: unphased.unwrap_or_default(),
        },
    )
    .parse(input)
}

/// Parses nucleotide alleles with uncertain phase.
///
/// Example: `123G>A(;)345del`
fn nucleotide_uncertain_allele(input: &str) -> ParseResult<'_, AlleleVariant<NucleotideEdit>> {
    let (input, first) = nucleotide_edit(input)?;
    let (input, _) = tag("(;)")(input)?;

    let (input, second) = alt((
        map(delimited(char('('), nucleotide_edit, char(')')), |edit| {
            Allele::uncertain_from_variants(vec![edit])
        }),
        map(nucleotide_edit, |edit| Allele::from_variants(vec![edit])),
    ))
    .parse(input)?;

    Ok((
        input,
        AlleleVariant {
            allele_one: Allele::from_variants(vec![first]),
            allele_two: Some(second),
            phase: Some(AllelePhase::Uncertain),
            variants_unphased: vec![],
        },
    ))
}

/// Parses nucleotide cis, trans and uncertain alleles.
///
/// This is the main entrance for handling possible nucleotide allele descriptions.
///
/// Examples: `[A;B]`, `[A];[B]`, `A(;)B`
pub(super) fn nucleotide_allele(input: &str) -> ParseResult<'_, AlleleVariant<NucleotideEdit>> {
    alt((
        // [123G>A];[345del](;)789dup
        nucleotide_trans_allele,
        // 123G>A(;)345del
        nucleotide_uncertain_allele,
        // [123G>A;345del]
        nucleotide_cis_allele,
    ))
    .parse(input)
}

/// Parses one nucleotide edit: location + edit.
///
/// Examples: `33038255C>A`, `4072_5145del`
pub(super) fn nucleotide_edit(input: &str) -> ParseResult<'_, NucleotideEdit> {
    map(
        pair(nucleotide_location, nucleotide_edit_kind),
        |(location, kind)| NucleotideEdit { location, kind },
    )
    .parse(input)
}

/// Parses a genomic (g.) description: an edit or an allele.
///
/// Examples: `g.33038255C>A`, `g.[123G>A;345del]`
pub(super) fn genomic_description(input: &str) -> ParseResult<'_, VariantDescription> {
    preceded(
        tag("g."),
        alt((
            // g.[123G>A;345del]
            map(nucleotide_allele, |variant| {
                VariantDescription::GenomicAllele(AlleleForm::Single(
                    variant.map_t(GenomicOutcome::from),
                ))
            }),
            // g.33038255C>A
            map(nucleotide_edit, |edit| {
                VariantDescription::Genomic(GenomicOutcome::from(edit))
            }),
        )),
    )
    .parse(input)
}

/// Parses various nucleotide edit families: substitution, deletion, etc.
///
/// Examples: `G>A`, `del`, `dup`, `insT`, `CAG[23]`
fn nucleotide_edit_kind(input: &str) -> ParseResult<'_, NucleotideEditKind> {
    alt((
        // =
        value(NucleotideEditKind::NoChange, char('=')),
        // delinsT, delinsN[12], delins[T;450_470;AGGG]
        map(
            preceded(tag("delins"), nucleotide_sequence_items),
            |items| NucleotideEditKind::DeletionInsertion { items },
        ),
        // del
        value(NucleotideEditKind::Deletion, tag("del")),
        // dup
        value(NucleotideEditKind::Duplication, tag("dup")),
        // insT, ins[T;450_470;AGGG]
        map(preceded(tag("ins"), nucleotide_sequence_items), |items| {
            NucleotideEditKind::Insertion { items }
        }),
        // inv
        value(NucleotideEditKind::Inversion, tag("inv")),
        // CAG[23], [14]. Keep this after insertion/deletion keywords.
        repeat_edits,
        // C>A
        map(
            pair(nucleotide_literal, preceded(char('>'), nucleotide_literal)),
            |(reference, alternate)| NucleotideEditKind::Substitution {
                reference,
                alternate,
            },
        ),
    ))
    .parse(input)
}

/// Parses inserted or replacement sequence items in an `ins` or `delins`
/// variant.
///
/// Examples: `T`, `[T;450_470;AGGG]`, `N[12]`
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
///
/// Examples: `AGGG`, `N[12]`, `450_470`, `NC_000022.10:g.35788169_35788352`
fn nucleotide_sequence_item(input: &str) -> ParseResult<'_, NucleotideSequenceItem> {
    alt((
        // N[12], CAG[23]
        map(
            alt((known_repeat_edit, unknown_repeat_edit)),
            NucleotideSequenceItem::Repeat,
        ),
        // 450_470, NC_000022.10:g.35788169_35788352
        map(sequence_segment, NucleotideSequenceItem::Copied),
        // T, AGGG
        map(nucleotide_literal, |value| {
            NucleotideSequenceItem::Literal(LiteralSequenceItem { value })
        }),
    ))
    .parse(input)
}

/// Parses a segment- or interval-type edit component that comes from either
/// local (current) or remote (other) reference source.
///
/// Examples: `850_900inv`, `NC_000022.10:g.35788169_35788352`
fn sequence_segment(input: &str) -> ParseResult<'_, CopiedSequenceItem> {
    alt((remote_sequence_segment, same_reference_sequence_segment)).parse(input)
}

/// Parses a current-reference segment such as `850_900inv`.
fn same_reference_sequence_segment(input: &str) -> ParseResult<'_, CopiedSequenceItem> {
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

/// Parses a sequence segment on a different (remote) reference.
///
/// Example: `NC_000022.10:g.35788169_35788352`
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

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;
    use nom::Parser;

    use super::*;
    use crate::model::{
        AlleleStateCertainty, CoordinateSystem, NucleotideEditKind, NucleotideSequenceItem,
        Quantity, RepeatEdit, RepeatSequenceUnit,
    };

    #[test]
    fn parses_nucleotide_edit_branches() {
        assert_eq!(
            all_consuming(nucleotide_edit_kind).parse("=").unwrap().1,
            NucleotideEditKind::NoChange
        );
        assert_eq!(
            all_consuming(nucleotide_edit_kind).parse("del").unwrap().1,
            NucleotideEditKind::Deletion
        );
        assert!(all_consuming(nucleotide_edit_kind).parse("delA").is_err());
        assert_eq!(
            all_consuming(nucleotide_edit_kind).parse("dup").unwrap().1,
            NucleotideEditKind::Duplication
        );
        assert_eq!(
            all_consuming(nucleotide_edit_kind).parse("inv").unwrap().1,
            NucleotideEditKind::Inversion
        );
        assert!(matches!(
            all_consuming(nucleotide_edit_kind).parse("C>A").unwrap().1,
            NucleotideEditKind::Substitution { .. }
        ));
        assert!(matches!(
            all_consuming(nucleotide_edit_kind).parse("insT").unwrap().1,
            NucleotideEditKind::Insertion { .. }
        ));
        assert!(matches!(
            all_consuming(nucleotide_edit_kind)
                .parse("delinsT")
                .unwrap()
                .1,
            NucleotideEditKind::DeletionInsertion { .. }
        ));
        assert!(matches!(
            all_consuming(nucleotide_edit_kind).parse("[4]").unwrap().1,
            NucleotideEditKind::Repeat { .. }
        ));
        assert!(matches!(
            all_consuming(nucleotide_edit_kind)
                .parse("CAG[23]")
                .unwrap()
                .1,
            NucleotideEditKind::Repeat { .. }
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
            NucleotideSequenceItem::Repeat(RepeatEdit {
                unit: Some(RepeatSequenceUnit::Unknown),
                quantity: Quantity::Known { count }
            }) if *count == 12
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
    fn parses_genomic_uncertain_allele_state() {
        let (_, variant) = all_consuming(nucleotide_allele)
            .parse("123G>A(;)345del")
            .unwrap();

        assert_eq!(variant.phase, Some(AllelePhase::Uncertain));

        assert_eq!(
            variant.allele_one.state_certainty,
            AlleleStateCertainty::Certain
        );

        assert_eq!(
            variant.allele_two.as_ref().unwrap().state_certainty,
            AlleleStateCertainty::Certain
        );

        let (_, variant) = all_consuming(nucleotide_allele)
            .parse("123G>A(;)(123G>A)")
            .unwrap();

        assert_eq!(variant.phase, Some(AllelePhase::Uncertain));

        assert_eq!(
            variant.allele_one.state_certainty,
            AlleleStateCertainty::Certain
        );

        assert_eq!(
            variant.allele_two.as_ref().unwrap().state_certainty,
            AlleleStateCertainty::Uncertain
        );
    }

    #[test]
    fn rejects_bracketed_genomic_uncertain_phase_forms() {
        assert!(all_consuming(nucleotide_allele)
            .parse("[123G>A](;)345del")
            .is_err());

        assert!(all_consuming(nucleotide_allele)
            .parse("[123G>A](;)(345del)")
            .is_err());
    }
}
