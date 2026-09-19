use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::char;
use nom::combinator::{map, opt, value};
use nom::multi::{many0, separated_list1};
use nom::sequence::{delimited, pair, preceded, separated_pair};
use nom::Parser;

use super::location::nucleotide_location;
use super::nucleotide::nucleotide_edit;
use super::repeat::{known_repeat_copy, known_repeat_unit, uncertain_repeat_copy};
use super::ParseResult;
use crate::model::{
    Allele, AlleleForm, AllelePhase, AlleleVariant, DerivedAllele, NucleotideEdit,
    NucleotideEditKind, OutcomeCertainty, Quantity, RepeatEdit, RnaOutcome, VariantDescription,
};

pub(super) fn produced_rna_outcome(input: &str) -> ParseResult<'_, RnaOutcome> {
    alt((
        // r.(A)
        map(delimited(char('('), nucleotide_edit, char(')')), |edit| {
            RnaOutcome::Produced {
                edit,
                certainty: OutcomeCertainty::Predicted,
            }
        }),
        // r.A
        map(nucleotide_edit, |edit| RnaOutcome::Produced {
            edit,
            certainty: OutcomeCertainty::Certain,
        }),
    ))
    .parse(input)
}

pub(super) fn special_rna_outcome(input: &str) -> ParseResult<'_, RnaOutcome> {
    alt((
        // `r.?`
        value(RnaOutcome::Unknown, char('?')),
        // `r.(?)`
        value(RnaOutcome::Indeterminate, tag("(?)")),
        // `r.0?`
        value(
            RnaOutcome::NoneProduced(OutcomeCertainty::Predicted),
            tag("0?"),
        ),
        // `r.0`
        value(
            RnaOutcome::NoneProduced(OutcomeCertainty::Certain),
            char('0'),
        ),
        // r.=
        value(RnaOutcome::NoChange(OutcomeCertainty::Certain), char('=')),
        // r.(=)
        value(
            RnaOutcome::NoChange(OutcomeCertainty::Predicted),
            tag("(=)"),
        ),
        // `r.spl?`, `r.spl`
        value(
            RnaOutcome::UncertainSplicing,
            alt((tag("spl?"), tag("spl"))),
        ),
    ))
    .parse(input)
}

pub(super) fn rna_outcome(input: &str) -> ParseResult<'_, RnaOutcome> {
    alt((special_rna_outcome, produced_rna_outcome)).parse(input)
}

pub(super) fn rna_variants_on_allele(input: &str) -> ParseResult<'_, Vec<RnaOutcome>> {
    alt((
        // (578c>u;1339a>g;1680del)
        map(
            delimited(
                char('('),
                separated_list1(char(';'), nucleotide_edit),
                char(')'),
            ),
            |edits| {
                edits
                    .into_iter()
                    .map(|edit| RnaOutcome::Produced {
                        edit,
                        certainty: OutcomeCertainty::Predicted,
                    })
                    .collect()
            },
        ),
        // 76a>u;103del
        // 76a>u;(103del)
        separated_list1(char(';'), produced_rna_outcome),
    ))
    .parse(input)
}

pub(super) fn rna_allele_component(input: &str) -> ParseResult<'_, Allele<RnaOutcome>> {
    map(
        delimited(
            char('['),
            alt((
                map(special_rna_outcome, |outcome| vec![outcome]),
                rna_variants_on_allele,
            )),
            char(']'),
        ),
        Allele::from_variants,
    )
    .parse(input)
}

pub(super) fn rna_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
    map(rna_allele_component, |allele| AlleleVariant {
        allele_one: allele,
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

pub(super) fn rna_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
    map(
        pair(
            separated_pair(rna_allele_component, char(';'), rna_allele_component),
            opt(preceded(
                tag("(;)"),
                separated_list1(tag("(;)"), rna_outcome),
            )),
        ),
        |((a1, a2), unphased)| AlleleVariant {
            allele_one: a1,
            allele_two: Some(a2),
            phase: Some(AllelePhase::Trans),
            variants_unphased: unphased.unwrap_or_default(),
        },
    )
    .parse(input)
}

pub(super) fn rna_uncertain_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
    let (input, a1) = produced_rna_outcome(input)?;
    let (input, _) = tag("(;)")(input)?;

    let (input, a2) = alt((
        // A(;)(B)
        //
        // Parentheses here mark uncertainty of the second allele state,
        // not prediction of the RNA outcome itself.
        map(delimited(char('('), nucleotide_edit, char(')')), |edit| {
            Allele::uncertain_from_variants(vec![RnaOutcome::Produced {
                edit,
                certainty: OutcomeCertainty::Certain,
            }])
        }),
        // A(;)B
        map(produced_rna_outcome, |outcome| {
            Allele::from_variants(vec![outcome])
        }),
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

// r.[location][14];[18], r.[position][unit][14];[18]
pub(super) fn rna_repeat_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
    let (input, location) = nucleotide_location(input)?;

    let (input, (unit, (q1, q2))) = pair(
        opt(known_repeat_unit),
        separated_pair(
            alt((uncertain_repeat_copy, known_repeat_copy)),
            char(';'),
            alt((uncertain_repeat_copy, known_repeat_copy)),
        ),
    )
    .parse(input)?;

    if unit.is_some() && !location.is_pos() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    let get_certainty = |q: &Quantity| match q {
        Quantity::Uncertain { .. } => OutcomeCertainty::Predicted,
        _ => OutcomeCertainty::Certain,
    };
    let rpt_one_certainty = get_certainty(&q1);
    let rpt_two_certainty = get_certainty(&q2);

    let rpt_one = RnaOutcome::Produced {
        edit: NucleotideEdit {
            location: location.clone(),
            kind: NucleotideEditKind::Repeat {
                blocks: vec![RepeatEdit {
                    quantity: q1,
                    unit: unit.clone(),
                }],
            },
        },
        certainty: rpt_one_certainty,
    };

    let rpt_two = RnaOutcome::Produced {
        edit: NucleotideEdit {
            location: location.clone(),
            kind: NucleotideEditKind::Repeat {
                blocks: vec![RepeatEdit {
                    quantity: q2,
                    unit: unit.clone(),
                }],
            },
        },
        certainty: rpt_two_certainty,
    };

    Ok((
        input,
        AlleleVariant {
            allele_one: Allele::from_variants(vec![rpt_one]),
            allele_two: Some(Allele::from_variants(vec![rpt_two])),
            phase: Some(AllelePhase::Trans),
            variants_unphased: vec![],
        },
    ))
}

pub(super) fn rna_derived_allele_form(input: &str) -> ParseResult<'_, DerivedAllele<RnaOutcome>> {
    let (input, _) = char('[')(input)?;

    let (input, first) = produced_rna_outcome(input)?;
    let (input, _) = char(',')(input)?;
    let (input, second) = produced_rna_outcome(input)?;
    let (input, rest) = many0(preceded(char(','), produced_rna_outcome)).parse(input)?;

    let (input, _) = char(']')(input)?;

    let mut outcomes = Vec::with_capacity(2 + rest.len());
    outcomes.push(first);
    outcomes.push(second);
    outcomes.extend(rest);

    Ok((input, DerivedAllele::from_outcomes(outcomes)))
}

pub(super) fn rna_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
    alt((
        rna_repeat_trans_allele,
        rna_trans_allele,
        rna_uncertain_allele,
        rna_cis_allele,
    ))
    .parse(input)
}

pub(super) fn rna_description(input: &str) -> ParseResult<'_, VariantDescription> {
    preceded(
        tag("r."),
        alt((
            map(rna_derived_allele_form, |derived| {
                VariantDescription::RnaAllele(AlleleForm::Derived(derived))
            }),
            map(rna_allele, |variant| {
                VariantDescription::RnaAllele(AlleleForm::Single(variant))
            }),
            map(rna_outcome, VariantDescription::Rna),
        )),
    )
    .parse(input)
}

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;
    use nom::Parser;

    use super::*;
    use crate::model::{AlleleForm, AlleleStateCertainty};

    #[test]
    fn parses_rna_uncertain_allele_state() {
        let (_, variant) = all_consuming(rna_allele).parse("76a>u(;)(76a>u)").unwrap();

        assert_eq!(variant.phase, Some(AllelePhase::Uncertain));

        let a1 = &variant.allele_one;
        let a2 = variant.allele_two.as_ref().unwrap();

        assert_eq!(a1.state_certainty, AlleleStateCertainty::Certain);
        assert_eq!(a2.state_certainty, AlleleStateCertainty::Uncertain);

        assert!(matches!(
            a1.variants.first().unwrap(),
            RnaOutcome::Produced {
                certainty: OutcomeCertainty::Certain,
                ..
            }
        ));

        assert!(matches!(
            a2.variants.first().unwrap(),
            RnaOutcome::Produced {
                certainty: OutcomeCertainty::Certain,
                ..
            }
        ));
    }

    #[test]
    fn parses_predicted_rna_outcome() {
        let (_, outcome) = all_consuming(rna_outcome).parse("(76a>u)").unwrap();

        assert!(matches!(
            outcome,
            RnaOutcome::Produced {
                certainty: OutcomeCertainty::Predicted,
                ..
            }
        ));
    }

    #[test]
    fn parses_rna_derived_allele_form() {
        let (_, description) = all_consuming(rna_description)
            .parse("r.[897u>g,832_960del,950a>g]")
            .unwrap();

        let VariantDescription::RnaAllele(AlleleForm::Derived(derived)) = description else {
            panic!("expected derived RNA allele form");
        };

        assert_eq!(derived.outcomes.len(), 3);

        assert!(derived.outcomes.iter().all(|outcome| {
            matches!(
                outcome,
                RnaOutcome::Produced {
                    certainty: OutcomeCertainty::Certain,
                    ..
                }
            )
        }));
    }

    #[test]
    fn rejects_single_form_as_rna_derived_allele_form() {
        assert!(all_consuming(rna_derived_allele_form)
            .parse("[897u>g]")
            .is_err());
    }
}
