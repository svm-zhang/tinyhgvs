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

pub(super) fn nucleotide_variants_on_allele(input: &str) -> ParseResult<'_, Vec<NucleotideEdit>> {
    delimited(
        char('['),
        separated_list1(char(';'), nucleotide_edit),
        char(']'),
    )
    .parse(input)
}

// [A;B]
// [(A;B)]
fn nucleotide_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<NucleotideEdit>> {
    map(nucleotide_variants_on_allele, |variants| AlleleVariant {
        allele_one: Allele::from_variants(variants),
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

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

fn nucleotide_uncertain_allele(input: &str) -> ParseResult<'_, AlleleVariant<NucleotideEdit>> {
    map(
        separated_pair(nucleotide_edit, tag("(;)"), nucleotide_edit),
        |(a1, a2)| AlleleVariant {
            allele_one: Allele::from_variants(vec![a1]),
            allele_two: Some(Allele::from_variants(vec![a2])),
            phase: Some(AllelePhase::Uncertain),
            variants_unphased: vec![],
        },
    )
    .parse(input)
}

pub(super) fn nucleotide_allele(input: &str) -> ParseResult<'_, AlleleVariant<NucleotideEdit>> {
    alt((
        nucleotide_trans_allele,
        nucleotide_uncertain_allele,
        nucleotide_cis_allele,
    ))
    .parse(input)
}

pub(super) fn nucleotide_edit(input: &str) -> ParseResult<'_, NucleotideEdit> {
    map(
        pair(nucleotide_location, nucleotide_edit_kind),
        |(location, kind)| NucleotideEdit { location, kind },
    )
    .parse(input)
}

pub(super) fn genomic_description(input: &str) -> ParseResult<'_, VariantDescription> {
    preceded(
        tag("g."),
        alt((
            map(nucleotide_allele, |variant| {
                VariantDescription::GenomicAllele(AlleleForm::Single(
                    variant.map_t(GenomicOutcome::from),
                ))
            }),
            map(nucleotide_edit, |edit| {
                VariantDescription::Genomic(GenomicOutcome::from(edit))
            }),
        )),
    )
    .parse(input)
}

/// Parses the currently supported nucleotide edit families.
fn nucleotide_edit_kind(input: &str) -> ParseResult<'_, NucleotideEditKind> {
    alt((
        value(NucleotideEditKind::NoChange, char('=')),
        map(
            preceded(tag("delins"), nucleotide_sequence_items),
            |items| NucleotideEditKind::DeletionInsertion { items },
        ),
        value(NucleotideEditKind::Deletion, tag("del")),
        value(NucleotideEditKind::Duplication, tag("dup")),
        map(preceded(tag("ins"), nucleotide_sequence_items), |items| {
            NucleotideEditKind::Insertion { items }
        }),
        value(NucleotideEditKind::Inversion, tag("inv")),
        // has to put repeat pattern behind insertion and deletion
        // insN[(100_120)] will be mistaken as repeat edit. The inserted
        // sequence item is a repeat but the repeat unit is not "insN"
        repeat_edits,
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
        // map(sequence_repeat, NucleotideSequenceItem::Repeat),
        // This creates one possible concern that it allows N[80], N[(80-100)],
        // and N[?] as one item of the insertion/delins edit items. However,
        // the HGVS standard does not say it is invalid syntax either.
        // - NC_000006.11:g.10791926_10791927ins[NC_000004.11:g.106370094_106370420;A[26]]
        // - NC_000006.11:g.10791926_10791927ins[NC_000004.11:g.106370094_106370420;N[26]]
        map(
            alt((known_repeat_edit, unknown_repeat_edit)),
            NucleotideSequenceItem::Repeat,
        ),
        map(sequence_segment, NucleotideSequenceItem::Copied),
        map(nucleotide_literal, |value| {
            NucleotideSequenceItem::Literal(LiteralSequenceItem { value })
        }),
    ))
    .parse(input)
}

/// Parses a segment- or interval-type edit component that comes from either
/// local (current) or remote (other) reference source.
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

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;
    use nom::Parser;

    use super::*;
    use crate::model::{
        CoordinateSystem, NucleotideEditKind, NucleotideSequenceItem, Quantity, RepeatEdit,
        RepeatSequenceUnit,
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
}
