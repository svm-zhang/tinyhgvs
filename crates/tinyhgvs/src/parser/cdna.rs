//! Coding-DNA description, outcome, and allele parsers.

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::char;
use nom::combinator::{map, value};
use nom::multi::separated_list1;
use nom::sequence::{delimited, preceded, separated_pair};
use nom::Parser;

use super::nucleotide::nucleotide_edit;
use super::ParseResult;

use crate::model::{
    Allele, AlleleForm, AllelePhase, AlleleVariant, CodingDnaOutcome, VariantDescription,
};

/// Parses one or more coding-DNA outcomes written on the same allele.
///
/// Examples: `2376G>C`, `2376G>C;2376=`
fn cdna_variants_on_allele(input: &str) -> ParseResult<'_, Vec<CodingDnaOutcome>> {
    separated_list1(char(';'), cdna_outcome).parse(input)
}

/// Parses one coding-DNA allele component.
///
/// Examples: `[?]`, `[2376G>C]`, `[2376G>C;2376=]`
fn cdna_allele_component(input: &str) -> ParseResult<'_, Allele<CodingDnaOutcome>> {
    map(
        delimited(
            char('['),
            alt((
                // [?]
                map(char('?'), |_| vec![CodingDnaOutcome::Unknown]),
                // [2376G>C;2376=]
                cdna_variants_on_allele,
            )),
            char(']'),
        ),
        Allele::from_variants,
    )
    .parse(input)
}

/// Parses in-cis coding-DNA allele.
///
/// Example: `[2376G>C;2376=]`
fn cdna_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<CodingDnaOutcome>> {
    map(cdna_allele_component, |allele| AlleleVariant {
        allele_one: allele,
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

/// Parses coding-DNA alleles in trans.
///
/// Example: `[2376G>C];[2376=]`
fn cdna_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<CodingDnaOutcome>> {
    map(
        separated_pair(cdna_allele_component, char(';'), cdna_allele_component),
        |(a1, a2)| AlleleVariant {
            allele_one: a1,
            allele_two: Some(a2),
            phase: Some(AllelePhase::Trans),
            variants_unphased: vec![],
        },
    )
    .parse(input)
}

/// Parses coding-DNA alleles with uncertain phase.
///
/// Examples: `76A>G(;)80del`, `76A>G(;)(80del)`
fn cdna_uncertain_allele(input: &str) -> ParseResult<'_, AlleleVariant<CodingDnaOutcome>> {
    let (input, a1) = cdna_outcome(input)?;
    let (input, _) = tag("(;)")(input)?;

    let (input, a2) = alt((
        // A(;)(B)
        map(delimited(char('('), cdna_outcome, char(')')), |outcome| {
            Allele::uncertain_from_variants(vec![outcome])
        }),
        // A(;)B
        map(cdna_outcome, |outcome| Allele::from_variants(vec![outcome])),
    ))
    .parse(input)?;

    Ok((
        input,
        AlleleVariant {
            allele_one: Allele::from_variants(vec![a1]),
            allele_two: Some(a2),
            phase: Some(AllelePhase::Uncertain),
            variants_unphased: vec![],
        },
    ))
}

/// Parses in-cis, in-trans, and uncertain coding-DNA allele descriptions.
///
/// This is the main entrance for handling possible cdna allele descriptions.
///
/// Examples: `[A;B]`, `[A];[B]`, `A(;)B`, `A(;)(B)`
fn cdna_allele(input: &str) -> ParseResult<'_, AlleleVariant<CodingDnaOutcome>> {
    alt((
        // [A];[B]
        cdna_trans_allele,
        // A(;)B, A(;)(B)
        cdna_uncertain_allele,
        // [A;B]
        cdna_cis_allele,
    ))
    .parse(input)
}

/// Parses one coding-DNA outcome.
///
/// Examples: `?`, `357+1G>A`, `4072_5145del`
fn cdna_outcome(input: &str) -> ParseResult<'_, CodingDnaOutcome> {
    alt((
        // ?
        value(CodingDnaOutcome::Unknown, char('?')),
        // 357+1G>A, 4072_5145del
        map(nucleotide_edit, CodingDnaOutcome::Known),
    ))
    .parse(input)
}

/// Parses a coding-DNA (c.) description.
///
/// Examples: `c.357+1G>A`, `c.[2376G>C];[2376=]`
pub(super) fn cdna_description(input: &str) -> ParseResult<'_, VariantDescription> {
    preceded(
        tag("c."),
        alt((
            // c.[2376G>C];[2376=]
            map(cdna_allele, |variant| {
                VariantDescription::CodingDnaAllele(AlleleForm::Single(variant))
            }),
            // c.357+1G>A, c.?
            map(cdna_outcome, VariantDescription::CodingDna),
        )),
    )
    .parse(input)
}

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;
    use nom::Parser;

    use super::*;
    use crate::model::AlleleStateCertainty;

    #[test]
    fn parses_cdna_uncertain_allele_state() {
        // A(;)B
        let (_, variant) = all_consuming(cdna_allele).parse("76A>G(;)80del").unwrap();

        assert_eq!(variant.phase, Some(AllelePhase::Uncertain));
        assert_eq!(
            variant.allele_one.state_certainty,
            AlleleStateCertainty::Certain
        );
        assert_eq!(
            variant.allele_two.as_ref().unwrap().state_certainty,
            AlleleStateCertainty::Certain
        );

        // A(;)(B)
        let (_, variant) = all_consuming(cdna_allele).parse("76A>G(;)(80del)").unwrap();

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
}
