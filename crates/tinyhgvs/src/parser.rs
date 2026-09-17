//! HGVS variant parsers.
//!
//! The parsing strategy:
//!
//! - parse the accepted syntax into the Rust data model
//! - route rejected inputs to the lightweight diagnostic classifier

use nom::branch::alt;
use nom::bytes::complete::{tag, take_while1};
use nom::character::complete::{char, digit1, one_of};
use nom::combinator::{all_consuming, map, map_res, opt, value, verify};
use nom::multi::{many1, separated_list1};
use nom::sequence::{delimited, pair, preceded, separated_pair, terminated};
use nom::{IResult, Parser};

use crate::diagnostics::classify_parse_failure;
use crate::error::ParseHgvsError;
use crate::model::{
    Accession, Allele, AllelePhase, AlleleVariant, CodingDnaOutcome, CoordinateSystem,
    CopiedSequenceItem, GenomicOutcome, HgvsVariant, Interval, LiteralSequenceItem, Location,
    NucleotideAnchor, NucleotideCoordinate, NucleotideEdit, NucleotideEditKind,
    NucleotideSequenceItem, OutcomeCertainty, ProteinCoordinate, ProteinEdit, ProteinEditKind,
    ProteinExtensionEdit, ProteinExtensionTerminal, ProteinFrameshiftStop,
    ProteinFrameshiftStopKind, ProteinOutcome, ProteinSequence, Quantity, ReferenceSpec,
    RepeatEdit, RepeatSequenceUnit, RnaOutcome, VariantDescription,
};

type ParseResult<'a, T> = IResult<&'a str, T>;

// Same numeric parser, named by semantic role at the frameshift call site.
use self::parse_quantity as parse_stop_ordinal;

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
///         assert_eq!(nucleotide.location.start().unwrap().anchor().unwrap(), NucleotideAnchor::Absolute);
///         assert_eq!(nucleotide.location.start().unwrap().coordinate().unwrap(), 357);
///         assert_eq!(nucleotide.location.start().unwrap().offset().unwrap(), 1);
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
///         assert_eq!(nucleotide.location.start().unwrap().anchor().unwrap(), NucleotideAnchor::RelativeCdsStart);
///         assert_eq!(nucleotide.location.start().unwrap().coordinate().unwrap(), -1);
///         assert_eq!(nucleotide.location.start().unwrap().offset().unwrap(), 0);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
///
/// CDS-anchored intronic positions in the 5' and 3' UTR:
///
/// ```rust
/// use tinyhgvs::{NucleotideAnchor, VariantDescription, parse_hgvs};
///
/// let five_prime_intronic = parse_hgvs("NM_001385026.1:c.-106+2T>A").unwrap();
/// let three_prime_intronic = parse_hgvs("NM_001272071.2:c.*639-1G>A").unwrap();
///
/// match five_prime_intronic.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         assert_eq!(nucleotide.location.start().unwrap().anchor().unwrap(), NucleotideAnchor::RelativeCdsStart);
///         assert_eq!(nucleotide.location.start().unwrap().coordinate().unwrap(), -106);
///         assert_eq!(nucleotide.location.start().unwrap().offset().unwrap(), 2);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
///
/// match three_prime_intronic.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         assert_eq!(nucleotide.location.start().unwrap().anchor().unwrap(), NucleotideAnchor::RelativeCdsEnd);
///         assert_eq!(nucleotide.location.start().unwrap().coordinate().unwrap(), 639);
///         assert_eq!(nucleotide.location.start().unwrap().offset().unwrap(), -1);
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
///         assert!(matches!(protein.effect, ProteinEffect::Known { .. }));
///     }
///     _ => unreachable!("expected protein variant"),
/// }
/// ```
///
/// A repeated sequence is returned as a repeat edit:
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, RepeatEdit, Quantity, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_004006.3:r.-124_-123[14]").unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         let NucleotideEdit::Repeat { blocks } = nucleotide.edit else {
///             unreachable!("expected repeat edit");
///         };
///         assert_eq!(blocks, &[RepeatEdit {
///             unit: None, quantity: Quantity::Known {count: 14}
///         }]);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
/// ```
///
/// A nucleotide allele variant with two in-trans alleles:
///
/// ```rust
/// use tinyhgvs::{AllelePhase, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_004006.2:c.[2376G>C];[2376=]").unwrap();
///
/// match variant.description {
///     VariantDescription::NucleotideAllele(allele) => {
///         assert_eq!(allele.allele_one.variants.len(), 1);
///         assert!(allele.allele_two.is_some());
///         assert_eq!(allele.phase, Some(AllelePhase::Trans));
///     }
///     _ => unreachable!("expected nucleotide allele"),
/// }
/// ```
///
/// A protein frameshift can be parsed in either short or long form:
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(protein) => match protein.effect {
///         ProteinEffect::Known { edit: ProteinEdit::Frameshift { to_residue, stop }, .. } => {
///             assert_eq!(to_residue.as_deref(), Some("Pro"));
///             assert_eq!(stop.ordinal, Some(23));
///         }
///         _ => unreachable!("expected protein frameshift"),
///     },
///     _ => unreachable!("expected protein variant"),
/// }
/// ```
///
/// A protein extension keeps the extended terminus, the first new residue when
/// present, and the new terminal ordinal together:
///
/// ```rust
/// use tinyhgvs::{ProteinEdit, ProteinEffect, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17").unwrap();
///
/// match variant.description {
///     VariantDescription::Protein(protein) => match protein.effect {
///         ProteinEffect::Known { edit: ProteinEdit::Extension(extension), .. } => {
///             assert_eq!(extension.to_residue.as_deref(), Some("Gln"));
///             assert_eq!(extension.terminal_ordinal, Some(17));
///         }
///         _ => unreachable!("expected protein extension"),
///     },
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
    alt((protein_variant, nucleotide_variant)).parse(input)
}

/// Parses the full HGVS variant (with a reference identifier).
fn nucleotide_variant(input: &str) -> ParseResult<'_, HgvsVariant> {
    // Parses the reference field
    let (input, reference) = terminated(reference_spec, char(':')).parse(input)?;
    let (input, (coordinate_system, description)) = alt((
        map(genomic_description, |description| {
            (CoordinateSystem::Genomic, description)
        }),
        map(cdna_description, |description| {
            (CoordinateSystem::CodingDna, description)
        }),
        map(rna_description, |description| {
            (CoordinateSystem::Rna, description)
        }),
    ))
    .parse(input)?;

    Ok((
        input,
        HgvsVariant::from(Some(reference), coordinate_system, description),
    ))
}

/// Parses context-dependent shorthand protein-level variant.
fn protein_variant(input: &str) -> ParseResult<'_, HgvsVariant> {
    // protein variant with sequence identifier field
    // intentionally leave the following as a local parser for easy following
    let with_sid = |i| map(terminated(reference_spec, char(':')), |reference| reference).parse(i);
    map(
        pair(terminated(opt(with_sid), tag("p.")), protein_description),
        |(refspec, description)| HgvsVariant::from(refspec, CoordinateSystem::Protein, description),
    )
    .parse(input)
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
        value(CoordinateSystem::CodingDna, char('c')),
        value(CoordinateSystem::Rna, char('r')),
        value(CoordinateSystem::Protein, char('p')),
    ))
    .parse(input)
}

// fn verify_failure<T>(input: &str) -> ParseResult<'_, T> {
//     Err(nom::Err::Error(nom::error::Error::new(
//         input,
//         nom::error::ErrorKind::Verify,
//     )))
// }

/// Parses a reusable range surface written as `thing_thing`.
fn range_with<T, P>(input: &str, parse_item: P) -> ParseResult<'_, Interval<T>>
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

/// Parser for phase marker written in allele variant description.
fn phase_marker(input: &str) -> ParseResult<'_, AllelePhase> {
    alt((
        // ;
        value(AllelePhase::Trans, char(';')),
        // (;)
        value(AllelePhase::Uncertain, tag("(;)")),
    ))
    .parse(input)
}

fn produced_rna_outcome(input: &str) -> ParseResult<'_, RnaOutcome> {
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

fn special_rna_outcome(input: &str) -> ParseResult<'_, RnaOutcome> {
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

fn rna_outcome(input: &str) -> ParseResult<'_, RnaOutcome> {
    alt((special_rna_outcome, produced_rna_outcome)).parse(input)
}

fn rna_variants_on_allele(input: &str) -> ParseResult<'_, Vec<RnaOutcome>> {
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

fn rna_allele_component(input: &str) -> ParseResult<'_, Allele<RnaOutcome>> {
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

fn rna_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
    map(rna_allele_component, |allele| AlleleVariant {
        allele_one: allele,
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

fn rna_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
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

fn rna_uncertain_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
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
fn rna_repeat_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
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

fn rna_allele(input: &str) -> ParseResult<'_, AlleleVariant<RnaOutcome>> {
    alt((
        rna_repeat_trans_allele,
        rna_trans_allele,
        rna_uncertain_allele,
        rna_cis_allele,
    ))
    .parse(input)
}

fn rna_description(input: &str) -> ParseResult<'_, VariantDescription> {
    preceded(
        tag("r."),
        alt((
            map(rna_allele, VariantDescription::RnaAllele),
            map(rna_outcome, VariantDescription::Rna),
        )),
    )
    .parse(input)
}

fn protein_outcome(input: &str) -> ParseResult<'_, ProteinOutcome> {
    alt((
        // (Ser68Arg)
        map(delimited(char('('), protein_edit, char(')')), |edit| {
            ProteinOutcome::Produced {
                edit,
                certainty: OutcomeCertainty::Predicted,
            }
        }),
        // Ser68Arg
        map(protein_edit, |edit| ProteinOutcome::Produced {
            edit,
            certainty: OutcomeCertainty::Certain,
        }),
    ))
    .parse(input)
}

fn special_protein_outcome(input: &str) -> ParseResult<'_, ProteinOutcome> {
    alt((
        // p.?
        value(ProteinOutcome::Unknown, char('?')),
        // p.0?
        value(
            ProteinOutcome::NoneProduced(OutcomeCertainty::Predicted),
            tag("0?"),
        ),
        // p.0
        value(
            ProteinOutcome::NoneProduced(OutcomeCertainty::Certain),
            char('0'),
        ),
    ))
    .parse(input)
}

fn protein_variants_on_allele(input: &str) -> ParseResult<'_, Vec<ProteinOutcome>> {
    alt((
        // (Ser68Arg;Asn594del)
        map(
            delimited(
                char('('),
                separated_list1(char(';'), protein_edit),
                char(')'),
            ),
            |edits| {
                edits
                    .into_iter()
                    .map(|edit| ProteinOutcome::Produced {
                        edit,
                        certainty: OutcomeCertainty::Predicted,
                    })
                    .collect()
            },
        ),
        // Ser68Arg;Asn594del
        // Phe233Leu;(Cys690Trp)
        separated_list1(char(';'), protein_outcome),
    ))
    .parse(input)
}

fn protein_allele_component(input: &str) -> ParseResult<'_, Allele<ProteinOutcome>> {
    map(
        delimited(
            char('['),
            alt((
                // [?]
                map(char('?'), |_| vec![ProteinOutcome::Unknown]),
                // [0]
                map(char('0'), |_| {
                    vec![ProteinOutcome::NoneProduced(OutcomeCertainty::Certain)]
                }),
                protein_variants_on_allele,
            )),
            char(']'),
        ),
        |variants| Allele::from_variants(variants),
    )
    .parse(input)
}

fn protein_cis_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    map(protein_allele_component, |allele| AlleleVariant {
        allele_one: allele,
        allele_two: None,
        phase: None,
        variants_unphased: vec![],
    })
    .parse(input)
}

fn protein_trans_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    map(
        separated_pair(
            protein_allele_component,
            char(';'),
            protein_allele_component,
        ),
        |(a1, a2)| AlleleVariant {
            allele_one: a1,
            allele_two: Some(a2),
            phase: Some(AllelePhase::Trans),
            variants_unphased: vec![],
        },
    )
    .parse(input)
}

fn protein_uncertain_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    map(
        separated_pair(protein_outcome, tag("(;)"), protein_outcome),
        |(a1, a2)| AlleleVariant {
            allele_one: Allele::from_variants(vec![a1]),
            allele_two: Some(Allele::from_variants(vec![a2])),
            phase: Some(AllelePhase::Uncertain),
            variants_unphased: vec![],
        },
    )
    .parse(input)
}

fn protein_allele(input: &str) -> ParseResult<'_, AlleleVariant<ProteinOutcome>> {
    alt((
        protein_trans_allele,
        protein_uncertain_allele,
        protein_cis_allele,
    ))
    .parse(input)
}

/// Parser for protein variant and allele description.
fn protein_description(input: &str) -> ParseResult<'_, VariantDescription> {
    alt((
        map(protein_allele, VariantDescription::ProteinAllele),
        map(special_protein_outcome, VariantDescription::Protein),
        map(protein_outcome, VariantDescription::Protein),
    ))
    .parse(input)
}

fn nucleotide_variants_on_allele(input: &str) -> ParseResult<'_, Vec<NucleotideEdit>> {
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
            // separated_pair(nucleotide_variants_on_allele, char(';'), |i| {
            //     allele_parser(i, nucleotide_variants_on_allele)
            // }),
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

fn nucleotide_allele(input: &str) -> ParseResult<'_, AlleleVariant<NucleotideEdit>> {
    alt((
        nucleotide_trans_allele,
        nucleotide_uncertain_allele,
        nucleotide_cis_allele,
    ))
    .parse(input)
}

fn genomic_description(input: &str) -> ParseResult<'_, VariantDescription> {
    preceded(
        tag("g."),
        alt((
            map(nucleotide_allele, |v| {
                VariantDescription::GenomicAllele(v.map_t(GenomicOutcome::from))
            }),
            map(nucleotide_edit, |v| {
                VariantDescription::Genomic(GenomicOutcome::Known(v))
            }),
        )),
    )
    .parse(input)
}

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

fn cdna_description(input: &str) -> ParseResult<'_, VariantDescription> {
    preceded(
        tag("c."),
        alt((
            map(cdna_allele, VariantDescription::CodingDnaAllele),
            map(cdna_outcome, VariantDescription::CodingDna),
        )),
    )
    .parse(input)
}

fn nucleotide_edit(input: &str) -> ParseResult<'_, NucleotideEdit> {
    map(
        pair(nucleotide_location, nucleotide_edit_kind),
        |(location, kind)| NucleotideEdit { location, kind },
    )
    .parse(input)
}

// fn known_cdna_outcome(input: &str) -> ParseResult<'_, CodingDnaOutcome> {
//     map(nucleotide_edit, CodingDnaOutcome::Known).parse(input)
// }
//
// fn unknown_cdna_outcome(input: &str) -> ParseResult<'_, CodingDnaOutcome> {
//     value(CodingDnaOutcome::Unknown, char('?')).parse(input)
// }

fn nucleotide_location(input: &str) -> ParseResult<'_, Location<NucleotideCoordinate>> {
    // Reject location description such as `(?_?)`, `(?_?)_(?_?)`.
    let is_valid_uncertain_location = |loc: &Interval<Interval<NucleotideCoordinate>>| {
        !(loc.start.is_fully_unknown() && loc.end.as_ref().map_or(true, Interval::is_fully_unknown))
    };
    // Reject unknown location such as `(?_B)`, `(A_?).
    let is_valid_known_interval = |interval: &Interval<NucleotideCoordinate>| {
        !(interval.has_unknown_bound() && !interval.is_fully_unknown())
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
/// - A
/// - A_B
/// - ?_B
/// - A_?
/// - ?_?
fn nucleotide_interval(input: &str) -> ParseResult<'_, Interval<NucleotideCoordinate>> {
    alt((
        // Interval coordinate, `A_B`
        |input| range_with(input, nucleotide_coordinate),
        // Single position coordinate, `A`
        map(nucleotide_coordinate, |start| Interval { start, end: None }),
    ))
    .parse(input)
}

/// Parses one uncertain interval unit (with parenthesis)
///
/// - (A_B)
/// - (A_?)
/// - (?_B)
fn nucleotide_uncertain_interval(input: &str) -> ParseResult<'_, Interval<NucleotideCoordinate>> {
    delimited(char('('), nucleotide_interval, char(')')).parse(input)
}

/// Parses nucleotide uncertain location. One location can be either one
/// or two uncertain interval units (separated by '_')
fn nucleotide_uncertain_location(
    input: &str,
) -> ParseResult<'_, Interval<Interval<NucleotideCoordinate>>> {
    alt((
        // (A_B)_(C_D)
        |input| range_with(input, nucleotide_uncertain_interval),
        // (A_B)
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
        // A, or A+offset, or A-offset, where A is the coordinate
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

/// Parses an unsigned decimal integer into `i32`.
fn parse_i32(input: &str) -> ParseResult<'_, i32> {
    map_res(digit1, str::parse::<i32>).parse(input)
}

fn parse_position(input: &str) -> ParseResult<'_, i32> {
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
fn parse_quantity(input: &str) -> ParseResult<'_, usize> {
    map_res(digit1, str::parse::<usize>).parse(input)
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

fn known_repeat_unit(input: &str) -> ParseResult<'_, RepeatSequenceUnit> {
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

fn unknown_repeat_unit(input: &str) -> ParseResult<'_, RepeatSequenceUnit> {
    value(RepeatSequenceUnit::Unknown, one_of("Nn")).parse(input)
}

fn known_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
    delimited(
        char('['),
        map(parse_quantity, |count| Quantity::Known { count }),
        char(']'),
    )
    .parse(input)
}

fn unknown_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
    delimited(char('['), value(Quantity::Unknown, char('?')), char(']')).parse(input)
}

fn uncertain_repeat_copy(input: &str) -> ParseResult<'_, Quantity> {
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

fn known_repeat_edit(input: &str) -> ParseResult<'_, RepeatEdit> {
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

fn unknown_repeat_edit(input: &str) -> ParseResult<'_, RepeatEdit> {
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
fn repeat_edits(input: &str) -> ParseResult<'_, NucleotideEditKind> {
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

/// Parses an edit component of literal nucleotide base changes.
fn nucleotide_literal(input: &str) -> ParseResult<'_, String> {
    map(
        take_while1(|c: char| c.is_ascii_alphabetic()),
        str::to_string,
    )
    .parse(input)
}

/// Parser for one known protein consequence, which means:
///
/// - `Ser68Arg`
/// - `Ser68del`
/// - `Ser68_Ala74insSerGln`
///
/// The parser does not digest `?` and `0`.
fn protein_edit(input: &str) -> ParseResult<'_, ProteinEdit> {
    map_res(
        pair(protein_location, protein_edit_kind),
        build_protein_edit_effect,
    )
    .parse(input)
}

fn build_protein_edit_effect(
    (location, kind): (Location<ProteinCoordinate>, ProteinEditKind),
) -> Result<ProteinEdit, ()> {
    let location = resolve_protein_effect_location(&location, &kind).ok_or(())?;
    Ok(ProteinEdit { location, kind })
}

fn resolve_protein_effect_location(
    location: &Location<ProteinCoordinate>,
    edit: &ProteinEditKind,
) -> Option<Location<ProteinCoordinate>> {
    let ProteinEditKind::Extension(extension) = edit else {
        return Some(location.clone());
    };

    let Location::Known(location) = location else {
        return None;
    };

    if location.end.is_some() {
        return None;
    }

    let mut start = location.start.clone();

    match extension.to_terminal {
        ProteinExtensionTerminal::N => {
            if start.residue != "Met"
                || start.ordinal != 1
                || extension.to_residue.is_some()
                || !matches!(extension.terminal_ordinal, Some(ordinal) if ordinal < 0)
            {
                return None;
            }
        }
        ProteinExtensionTerminal::C => {
            if start.residue != "Ter"
                || extension.to_residue.is_none()
                || matches!(extension.terminal_ordinal, Some(ordinal) if ordinal <= 0)
            {
                return None;
            }
            start.residue = "Ter".to_string();
        }
    }

    Some(Location::from_known(Interval { start, end: None }))
}

/// Parses the currently supported protein edit families.
fn protein_edit_kind(input: &str) -> ParseResult<'_, ProteinEditKind> {
    alt((
        value(
            ProteinEditKind::NoChange(OutcomeCertainty::Predicted),
            tag("(=)"),
        ),
        value(
            ProteinEditKind::NoChange(OutcomeCertainty::Certain),
            char('='),
        ),
        map(preceded(tag("delins"), protein_sequence), |sequence| {
            ProteinEditKind::DeletionInsertion { sequence }
        }),
        value(ProteinEditKind::Deletion, tag("del")),
        value(ProteinEditKind::Duplication, tag("dup")),
        protein_repeat,
        protein_extension_edit,
        protein_frameshift_edit,
        map(preceded(tag("ins"), protein_sequence), |sequence| {
            ProteinEditKind::Insertion { sequence }
        }),
        map(protein_symbol, |to| ProteinEditKind::Substitution { to }),
    ))
    .parse(input)
}

/// Parses one supported protein location, known or uncertain.
fn protein_location(input: &str) -> ParseResult<'_, Location<ProteinCoordinate>> {
    alt((
        // (Ala123_Pro131) and (Ala123_Pro131)_(Gly140_Leu142)
        map(protein_uncertain_location, Location::from_uncertain),
        // Trp24 and Lys23_Val25
        map(protein_interval, Location::from_known),
    ))
    .parse(input)
}

/**
Parses a single protein position or an interval.

- Single position: `Ala237`
- Interval: `Ala237_Pro161`
*/
fn protein_interval(input: &str) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    alt((
        |input| range_with(input, protein_coordinate),
        map(protein_coordinate, |start| Interval { start, end: None }),
    ))
    .parse(input)
}

/**
Parse one protein uncertain interval unit (with parenthesis). This is a wrapper
over protein_interval parser.

- `(Ala237_Pro161)`
*/
fn protein_uncertain_interval(input: &str) -> ParseResult<'_, Interval<ProteinCoordinate>> {
    delimited(char('('), protein_interval, char(')')).parse(input)
}

/**
Parses protein locations written with uncertain-region syntax.
*/
fn protein_uncertain_location(
    input: &str,
) -> ParseResult<'_, Interval<Interval<ProteinCoordinate>>> {
    alt((
        |input| range_with(input, protein_uncertain_interval),
        map(protein_uncertain_interval, |start| Interval {
            start,
            end: None,
        }),
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

fn protein_repeat(input: &str) -> ParseResult<'_, ProteinEditKind> {
    // p.Ala2[10]
    // p.(Gln18)[(70_80)]
    map(
        alt((known_repeat_copy, uncertain_repeat_copy)),
        |quantity| {
            ProteinEditKind::Repeat(RepeatEdit {
                unit: None,
                quantity,
            })
        },
    )
    .parse(input)
}

/// Parses N-terminal and C-terminal protein extension syntax.
fn protein_extension_edit(input: &str) -> ParseResult<'_, ProteinEditKind> {
    alt((
        map(
            preceded(tag("ext"), protein_n_terminal_extension_ordinal),
            |terminal_ordinal| {
                ProteinEditKind::Extension(ProteinExtensionEdit {
                    to_terminal: ProteinExtensionTerminal::N,
                    to_residue: None,
                    terminal_ordinal: Some(terminal_ordinal),
                })
            },
        ),
        map(
            pair(
                protein_extension_residue,
                protein_c_terminal_extension_state,
            ),
            |(to_residue, terminal_ordinal)| {
                ProteinEditKind::Extension(ProteinExtensionEdit {
                    to_terminal: ProteinExtensionTerminal::C,
                    to_residue: Some(to_residue),
                    terminal_ordinal,
                })
            },
        ),
    ))
    .parse(input)
}

/// Parses the required negative ordinal in N-terminal extension syntax.
fn protein_n_terminal_extension_ordinal(input: &str) -> ParseResult<'_, i32> {
    map(preceded(char('-'), parse_i32), |ordinal| -ordinal).parse(input)
}

/// Parses the residue replacing the reference stop codon in C-terminal extension syntax.
fn protein_extension_residue(input: &str) -> ParseResult<'_, String> {
    let (input, residue) = protein_symbol(input)?;

    if residue == "Ter" {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )))
    } else {
        Ok((input, residue))
    }
}

/// Parses the terminal state in C-terminal extension syntax.
fn protein_c_terminal_extension_state(input: &str) -> ParseResult<'_, Option<i32>> {
    preceded(
        tag("ext"),
        alt((
            value(None, pair(alt((tag("Ter"), tag("*"))), char('?'))),
            map(preceded(alt((tag("Ter"), tag("*"))), parse_i32), Some),
        )),
    )
    .parse(input)
}

/// Parses short and long protein frameshift syntax.
fn protein_frameshift_edit(input: &str) -> ParseResult<'_, ProteinEditKind> {
    alt((
        map(
            pair(
                protein_frameshift_residue,
                pair(tag("fs"), protein_frameshift_stop),
            ),
            |(to_residue, (_, stop))| ProteinEditKind::Frameshift {
                to_residue: Some(to_residue),
                stop,
            },
        ),
        value(
            ProteinEditKind::Frameshift {
                to_residue: None,
                stop: ProteinFrameshiftStop {
                    ordinal: None,
                    kind: ProteinFrameshiftStopKind::Omitted,
                },
            },
            tag("fs"),
        ),
    ))
    .parse(input)
}

/// Parses the explicit stop-state in long protein frameshift notation.
fn protein_frameshift_stop(input: &str) -> ParseResult<'_, ProteinFrameshiftStop> {
    alt((
        value(
            ProteinFrameshiftStop {
                ordinal: None,
                kind: ProteinFrameshiftStopKind::Unknown,
            },
            pair(alt((tag("Ter"), tag("*"))), char('?')),
        ),
        map(
            preceded(alt((tag("Ter"), tag("*"))), parse_stop_ordinal),
            |ordinal| ProteinFrameshiftStop {
                ordinal: Some(ordinal),
                kind: ProteinFrameshiftStopKind::Known,
            },
        ),
    ))
    .parse(input)
}

/// Parses the explicitly written first residue in long protein frameshift syntax.
fn protein_frameshift_residue(input: &str) -> ParseResult<'_, String> {
    let (input, residue) = protein_symbol(input)?;

    if residue == "Ter" {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )))
    } else {
        Ok((input, residue))
    }
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
            return Ok((rest, normalize_protein_symbol(symbol)));
        }
    }

    Err(nom::Err::Error(nom::error::Error::new(
        input,
        nom::error::ErrorKind::Tag,
    )))
}

fn normalize_protein_symbol(symbol: &str) -> String {
    if symbol == "*" {
        "Ter".to_string()
    } else {
        symbol.to_string()
    }
}

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;

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
    fn parses_special_protein_outcome_branches() {
        assert_eq!(
            all_consuming(special_protein_outcome).parse("?").unwrap().1,
            ProteinOutcome::Unknown
        );
        assert_eq!(
            all_consuming(special_protein_outcome).parse("0").unwrap().1,
            ProteinOutcome::NoneProduced(OutcomeCertainty::Certain)
        );
        assert_eq!(
            all_consuming(special_protein_outcome).parse("0?").unwrap().1,
            ProteinOutcome::NoneProduced(OutcomeCertainty::Predicted)
        );
    }

    #[test]
    fn parses_protein_outcome_branches() {
        assert!(matches!(
            all_consuming(protein_outcome).parse("(Trp24Ter)").unwrap().1,
            ProteinOutcome::Produced {
                edit: ProteinEdit {
                    kind: ProteinEditKind::Substitution { .. },
                    ..
                },
                certainty: OutcomeCertainty::Predicted,
            }
        ));
        assert!(matches!(
            all_consuming(protein_outcome).parse("Trp24Ter").unwrap().1,
            ProteinOutcome::Produced {
                edit: ProteinEdit {
                    kind: ProteinEditKind::Substitution { .. },
                    ..
                },
                certainty: OutcomeCertainty::Certain,
            }
        ));
        assert!(matches!(
            all_consuming(protein_outcome).parse("Ala2[10]").unwrap().1,
            ProteinOutcome::Produced {
                edit: ProteinEdit {
                    kind: ProteinEditKind::Repeat(RepeatEdit {
                        unit: None,
                        quantity: Quantity::Known { count: 10 },
                    }),
                    ..
                },
                certainty: OutcomeCertainty::Certain,
            }
        ));
        assert!(matches!(
            all_consuming(protein_outcome).parse("Arg97fs").unwrap().1,
            ProteinOutcome::Produced {
                edit: ProteinEdit {
                    kind: ProteinEditKind::Frameshift {
                        to_residue: None,
                        stop: ProteinFrameshiftStop {
                            ordinal: None,
                            kind: ProteinFrameshiftStopKind::Omitted,
                        },
                    },
                    ..
                },
                certainty: OutcomeCertainty::Certain,
            }
        ));
        assert!(matches!(
            all_consuming(protein_outcome).parse("Met1ext-5").unwrap().1,
            ProteinOutcome::Produced {
                edit: ProteinEdit {
                    kind: ProteinEditKind::Extension(ProteinExtensionEdit {
                        to_terminal: ProteinExtensionTerminal::N,
                        to_residue: None,
                        terminal_ordinal: Some(-5),
                    }),
                    ..
                },
                certainty: OutcomeCertainty::Certain,
            }
        ));
        assert!(matches!(
            all_consuming(protein_outcome)
                .parse("Ter110GlnextTer17")
                .unwrap()
                .1,
            ProteinOutcome::Produced {
                edit: ProteinEdit {
                    kind: ProteinEditKind::Extension(ProteinExtensionEdit {
                        to_terminal: ProteinExtensionTerminal::C,
                        to_residue: Some(_),
                        terminal_ordinal: Some(17),
                    }),
                    ..
                },
                certainty: OutcomeCertainty::Certain,
            }
        ));
        assert!(matches!(
            all_consuming(protein_outcome)
                .parse("Arg97ProfsTer23")
                .unwrap()
                .1,
            ProteinOutcome::Produced {
                edit: ProteinEdit {
                    kind: ProteinEditKind::Frameshift {
                        to_residue: Some(_),
                        stop: ProteinFrameshiftStop {
                            ordinal: Some(23),
                            kind: ProteinFrameshiftStopKind::Known,
                        },
                    },
                    ..
                },
                certainty: OutcomeCertainty::Certain,
            }
        ));
        assert!(matches!(
            all_consuming(protein_outcome)
                .parse("Ile327Argfs*?")
                .unwrap()
                .1,
            ProteinOutcome::Produced {
                edit: ProteinEdit {
                    kind: ProteinEditKind::Frameshift {
                        to_residue: Some(_),
                        stop: ProteinFrameshiftStop {
                            ordinal: None,
                            kind: ProteinFrameshiftStopKind::Unknown,
                        },
                    },
                    ..
                },
                certainty: OutcomeCertainty::Certain,
            }
        ));
    }

    #[test]
    fn parses_protein_frameshift_branches() {
        assert_eq!(
            all_consuming(protein_edit_kind).parse("fs").unwrap().1,
            ProteinEditKind::Frameshift {
                to_residue: None,
                stop: ProteinFrameshiftStop {
                    ordinal: None,
                    kind: ProteinFrameshiftStopKind::Omitted,
                },
            }
        );
        assert_eq!(
            all_consuming(protein_edit_kind).parse("ProfsTer23").unwrap().1,
            ProteinEditKind::Frameshift {
                to_residue: Some("Pro".to_string()),
                stop: ProteinFrameshiftStop {
                    ordinal: Some(23),
                    kind: ProteinFrameshiftStopKind::Known,
                },
            }
        );
        assert_eq!(
            all_consuming(protein_edit_kind).parse("Argfs*?").unwrap().1,
            ProteinEditKind::Frameshift {
                to_residue: Some("Arg".to_string()),
                stop: ProteinFrameshiftStop {
                    ordinal: None,
                    kind: ProteinFrameshiftStopKind::Unknown,
                },
            }
        );
        assert!(all_consuming(protein_edit_kind).parse("TerfsTer2").is_err());
    }

    #[test]
    fn parses_protein_extension_branches() {
        assert_eq!(
            all_consuming(protein_edit_kind).parse("ext-5").unwrap().1,
            ProteinEditKind::Extension(ProteinExtensionEdit {
                to_terminal: ProteinExtensionTerminal::N,
                to_residue: None,
                terminal_ordinal: Some(-5),
            })
        );
        assert_eq!(
            all_consuming(protein_edit_kind).parse("GlnextTer17").unwrap().1,
            ProteinEditKind::Extension(ProteinExtensionEdit {
                to_terminal: ProteinExtensionTerminal::C,
                to_residue: Some("Gln".to_string()),
                terminal_ordinal: Some(17),
            })
        );
        assert_eq!(
            all_consuming(protein_edit_kind).parse("Argext*?").unwrap().1,
            ProteinEditKind::Extension(ProteinExtensionEdit {
                to_terminal: ProteinExtensionTerminal::C,
                to_residue: Some("Arg".to_string()),
                terminal_ordinal: None,
            })
        );
        assert!(all_consuming(protein_edit_kind)
            .parse("TerextTer17")
            .is_err());
    }
}
