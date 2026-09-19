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

fn cdna_variants_on_allele(input: &str) -> ParseResult<'_, Vec<CodingDnaOutcome>> {
    separated_list1(char(';'), cdna_outcome).parse(input)
}

fn cdna_allele_component(input: &str) -> ParseResult<'_, Allele<CodingDnaOutcome>> {
    map(
        delimited(
            char('['),
            alt((
                map(char('?'), |_| vec![CodingDnaOutcome::Unknown]),
                cdna_variants_on_allele,
            )),
            char(']'),
        ),
        Allele::from_variants,
    )
    .parse(input)
}

fn cdna_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<CodingDnaOutcome>> {
    map(cdna_allele_component, |allele| AlleleVariant {
        allele_one: allele,
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

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

fn cdna_allele(input: &str) -> ParseResult<'_, AlleleVariant<CodingDnaOutcome>> {
    alt((cdna_trans_allele, cdna_uncertain_allele, cdna_cis_allele)).parse(input)
}

fn cdna_outcome(input: &str) -> ParseResult<'_, CodingDnaOutcome> {
    alt((
        value(CodingDnaOutcome::Unknown, char('?')),
        map(nucleotide_edit, CodingDnaOutcome::Known),
    ))
    .parse(input)
}

pub(super) fn cdna_description(input: &str) -> ParseResult<'_, VariantDescription> {
    preceded(
        tag("c."),
        alt((
            map(cdna_allele, |variant| {
                VariantDescription::CodingDnaAllele(AlleleForm::Single(variant))
            }),
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
