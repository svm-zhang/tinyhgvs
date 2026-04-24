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
use nom::multi::{many0, many1, separated_list1};
use nom::sequence::{delimited, pair, preceded, separated_pair};
use nom::{IResult, Parser};

use crate::diagnostics::classify_parse_failure;
use crate::error::ParseHgvsError;
use crate::model::{
    Accession, Allele, AllelePhase, AlleleVariant, CoordinateKind, CoordinateSystem,
    CopiedSequenceItem, HgvsVariant, Interval, LiteralSequenceItem, Location, NucleotideAnchor,
    NucleotideCoordinate, NucleotideEdit, NucleotideRepeatBlock, NucleotideSequenceItem,
    NucleotideVariant, ProteinCoordinate, ProteinEdit, ProteinEffect, ProteinExtensionEdit,
    ProteinExtensionTerminal, ProteinFrameshiftStop, ProteinFrameshiftStopKind, ProteinSequence,
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
///         assert_eq!(nucleotide.location.start().unwrap().anchor, NucleotideAnchor::Absolute);
///         assert_eq!(nucleotide.location.start().unwrap().coordinate, 357);
///         assert_eq!(nucleotide.location.start().unwrap().offset, 1);
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
///         assert_eq!(nucleotide.location.start().unwrap().anchor, NucleotideAnchor::RelativeCdsStart);
///         assert_eq!(nucleotide.location.start().unwrap().coordinate, -1);
///         assert_eq!(nucleotide.location.start().unwrap().offset, 0);
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
///         assert_eq!(nucleotide.location.start().unwrap().anchor, NucleotideAnchor::RelativeCdsStart);
///         assert_eq!(nucleotide.location.start().unwrap().coordinate, -106);
///         assert_eq!(nucleotide.location.start().unwrap().offset, 2);
///     }
///     _ => unreachable!("expected nucleotide variant"),
/// }
///
/// match three_prime_intronic.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         assert_eq!(nucleotide.location.start().unwrap().anchor, NucleotideAnchor::RelativeCdsEnd);
///         assert_eq!(nucleotide.location.start().unwrap().coordinate, 639);
///         assert_eq!(nucleotide.location.start().unwrap().offset, -1);
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
/// A repeated sequence is returned as a repeat edit:
///
/// ```rust
/// use tinyhgvs::{NucleotideEdit, VariantDescription, parse_hgvs};
///
/// let variant = parse_hgvs("NM_004006.3:r.-124_-123[14]").unwrap();
///
/// match variant.description {
///     VariantDescription::Nucleotide(nucleotide) => {
///         let NucleotideEdit::Repeat { blocks } = nucleotide.edit else {
///             unreachable!("expected repeat edit");
///         };
///         assert_eq!(blocks[0].count, 14);
///         assert_eq!(blocks[0].unit, None);
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
///         ProteinEffect::Edit { edit: ProteinEdit::Frameshift { to_residue, stop }, .. } => {
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
///         ProteinEffect::Edit { edit: ProteinEdit::Extension(extension), .. } => {
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
        nucleotide_description(coordinate_system, input)
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

/// Parser for nucleotide variant and allele description.
fn nucleotide_description(
    coordinate_system: CoordinateSystem,
    input: &str,
) -> ParseResult<'_, VariantDescription> {
    alt((
        // Nucleotide allele description: [123G>A;345del]
        // Composed of nucleotide variant descriptions
        |input| {
            allele_description(
                input,
                |input| nucleotide_initial_allele(coordinate_system, input),
                |phase, input| next_nucleotide_allele(coordinate_system, phase, input),
                VariantDescription::NucleotideAllele,
            )
        },
        // Nucleotide variant description: 123G>A
        map(
            |input| nucleotide_variant_description(coordinate_system, input),
            VariantDescription::Nucleotide,
        ),
    ))
    .parse(input)
}

fn verify_failure<T>(input: &str) -> ParseResult<'_, T> {
    Err(nom::Err::Error(nom::error::Error::new(
        input,
        nom::error::ErrorKind::Verify,
    )))
}

/// Universal parser for nucleotide and protein allele description carrying
/// variants, such as:
/// - `NC_000001.11:g.[123G>A;345del]`
/// - `NM_004006.3:r.[123c>a;345del]`
/// - `p.[Ser73Arg;Asn103del]`
/// - `NP_003997.1:p.[Ser68Arg;Asn594del]`
/// - `p.[(Ser73Arg;Asn103del)]`
fn variants_on_allele<T, PVariant>(input: &str, parse_variant: PVariant) -> ParseResult<'_, Vec<T>>
where
    PVariant: Fn(&str) -> ParseResult<'_, T>,
{
    separated_list1(char(';'), |input| parse_variant(input)).parse(input)
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

/// Main parser entry point for parsing both nucleotide and protein allele
/// description. The work sequence flow is:
/// - initial allele
/// - optional second established allele
/// - additional unphased alleles after `(;)`
fn allele_description<T, PInitial, PNext, BuildDescription>(
    input: &str,
    parse_initial_allele: PInitial,
    parse_next_allele: PNext,
    build_description: BuildDescription,
) -> ParseResult<'_, VariantDescription>
where
    PInitial: Fn(&str) -> ParseResult<'_, (Allele<T>, bool)>,
    PNext: Fn(AllelePhase, &str) -> ParseResult<'_, Allele<T>>,
    BuildDescription: Fn(AlleleVariant<T>) -> VariantDescription,
{
    // Initial allele
    let (input, (allele_one, initial_bracketed)) = parse_initial_allele(input)?;
    let allele_one_variant_count = allele_one.variants.len();

    // Optional second allele when present
    let (input, established_second) = opt(|input| {
        let (input, phase) = phase_marker(input)?;

        // Reject A;[B]
        if phase == AllelePhase::Trans && !initial_bracketed {
            return verify_failure(input);
        }

        // Reject [A](;)B
        if phase == AllelePhase::Uncertain && initial_bracketed && allele_one_variant_count == 1 {
            return verify_failure(input);
        }

        let (input, allele_two) = parse_next_allele(phase, input)?;
        Ok((input, (phase, allele_two)))
    })
    .parse(input)?;

    // Additional alleles
    let (input, alleles_unphased) = if established_second.is_some() {
        many0(|input| {
            let (input, _) = tag("(;)")(input)?;
            parse_next_allele(AllelePhase::Uncertain, input)
        })
        .parse(input)?
    } else {
        (input, Vec::new())
    };

    let (allele_two, phase) = match established_second {
        Some((phase, allele_two)) => (Some(allele_two), Some(phase)),
        None => (None, None),
    };

    // Reject A or A(;) without a written second allele.
    if !initial_bracketed && allele_two.is_none() {
        return verify_failure(input);
    }

    Ok((
        input,
        build_description(AlleleVariant {
            allele_one,
            allele_two,
            phase,
            alleles_unphased,
        }),
    ))
}

/// Parses one supported nucleotide variant description such as `123G>A` or
/// `357+1G>A`, without wrapping it in the top-level description enum.
fn nucleotide_variant_description(
    coordinate_system: CoordinateSystem,
    input: &str,
) -> ParseResult<'_, NucleotideVariant> {
    let (input, initial_location) = nucleotide_location(coordinate_system, input)?;
    let (input, edit) = nucleotide_edit(input)?;

    if !is_valid_nucleotide_repeat(coordinate_system, &initial_location, &edit) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    let Some(location) = resolve_nucleotide_location(&initial_location, &edit) else {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    };

    Ok((input, NucleotideVariant { location, edit }))
}

/// Parses one supported nucleotide location, known or uncertain.
fn nucleotide_location(
    coordinate_system: CoordinateSystem,
    input: &str,
) -> ParseResult<'_, Location<NucleotideCoordinate>> {
    match coordinate_system {
        CoordinateSystem::Rna => {
            map_res(nucleotide_interval, build_known_nucleotide_location).parse(input)
        }
        CoordinateSystem::Genomic
        | CoordinateSystem::CircularGenomic
        | CoordinateSystem::Mitochondrial
        | CoordinateSystem::CodingDna
        | CoordinateSystem::NonCodingDna => alt((
            // (71_72) and (123_234)_(345_456), (?_87), (123_?)_(?_456)
            map_res(
                nucleotide_uncertain_location,
                build_uncertain_nucleotide_location,
            ),
            // 93 and 93_94
            map_res(nucleotide_interval, build_known_nucleotide_location),
        ))
        .parse(input),
        CoordinateSystem::Protein => {
            unreachable!("nucleotide location on protein coordinate system")
        }
    }
}

fn build_known_nucleotide_location(
    interval: Interval<NucleotideCoordinate>,
) -> Result<Location<NucleotideCoordinate>, ()> {
    // Reject known-location branches that still contain `?` on either side.
    if interval.has_unknown_bound() {
        return Err(());
    }

    // Build one known nucleotide location from a fully known interval.
    Ok(Location::from_known(interval))
}

fn build_uncertain_nucleotide_location(
    location: Interval<Interval<NucleotideCoordinate>>,
) -> Result<Location<NucleotideCoordinate>, ()> {
    // Reject fully unknown uncertain regions such as `(?_?)` and `(?_?)_(A_B)`.
    if location.start.is_fully_unknown()
        || location
            .end
            .as_ref()
            .map_or(false, Interval::is_fully_unknown)
    {
        return Err(());
    }

    // Build one uncertain nucleotide location from parenthesized regions.
    Ok(Location::from_uncertain(location))
}

/// Parser for digesting the initial allele written in a nucleotide allele
/// description, such as
///
/// - `c.[123G>A;345del]`
/// - `c.123G>A(;)345del`
///
/// Returns both the parsed allele and whether that first allele was bracketed.
fn nucleotide_initial_allele(
    coordinate_system: CoordinateSystem,
    input: &str,
) -> ParseResult<'_, (Allele<NucleotideVariant>, bool)> {
    alt((
        // [123G>A;345del]
        map(
            |input| nucleotide_bracket_allele(coordinate_system, input),
            |allele| (allele, true),
        ),
        // 123G>A
        map(
            |input| nucleotide_variant_description(coordinate_system, input),
            |variant| (Allele::from_variants(vec![variant]), false),
        ),
    ))
    .parse(input)
}

/// Parser for pattern of nucleotide allele wrapped in bracket, such as:
///
/// - `[123G>A;345del]`
fn nucleotide_bracket_allele(
    coordinate_system: CoordinateSystem,
    input: &str,
) -> ParseResult<'_, Allele<NucleotideVariant>> {
    map(
        delimited(
            char('['),
            |input| {
                variants_on_allele(input, |input| {
                    nucleotide_variant_description(coordinate_system, input)
                })
            },
            char(']'),
        ),
        Allele::from_variants,
    )
    .parse(input)
}

/// Parser for digesting the next nucleotide allele following the initial
/// allele + phase marker in the nucleotide allele description, such as:
///
/// - `[345del]` in `NC_000001.11:g.[123G>A];[345del]` after `;`
/// - `3103del` in `NM_004006.2:c.2376G>C(;)3103del` after '(;)'
///
/// The phase marker after the initial allele imposes syntactic rule.
///
/// Examples:
/// - `NM_004006.2:c.2376G>C;3103del`: invalid
/// - `NM_004006.2:c.[2376G>C](;)[3103del]`: invalid
fn next_nucleotide_allele(
    coordinate_system: CoordinateSystem,
    phase: AllelePhase,
    input: &str,
) -> ParseResult<'_, Allele<NucleotideVariant>> {
    match phase {
        AllelePhase::Trans => nucleotide_bracket_allele(coordinate_system, input),
        AllelePhase::Uncertain => map(
            |input| nucleotide_variant_description(coordinate_system, input),
            |variant| Allele::from_variants(vec![variant]),
        )
        .parse(input),
    }
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
        map(
            separated_pair(nucleotide_coordinate, char('_'), nucleotide_coordinate),
            |(start, end)| Interval {
                start,
                end: Some(end),
            },
        ),
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
        map(
            separated_pair(
                nucleotide_uncertain_interval,
                char('_'),
                nucleotide_uncertain_interval,
            ),
            |(start, end)| Interval {
                start,
                end: Some(end),
            },
        ),
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
                preceded(char('-'), parse_nonzero_i32),
                opt(pair(alt((char('+'), char('-'))), parse_i32)),
            ),
            |(coordinate, offset)| {
                let offset = offset
                    .map(|(sign, value)| if sign == '-' { -value } else { value })
                    .unwrap_or(0);

                NucleotideCoordinate {
                    anchor: NucleotideAnchor::RelativeCdsStart,
                    coordinate: CoordinateKind::Known(-coordinate),
                    offset,
                }
            },
        ),
        map(
            pair(
                preceded(char('*'), parse_nonzero_i32),
                opt(pair(alt((char('+'), char('-'))), parse_i32)),
            ),
            |(coordinate, offset)| {
                let offset = offset
                    .map(|(sign, value)| if sign == '-' { -value } else { value })
                    .unwrap_or(0);

                NucleotideCoordinate {
                    anchor: NucleotideAnchor::RelativeCdsEnd,
                    coordinate: CoordinateKind::Known(coordinate),
                    offset,
                }
            },
        ),
        // A, or A+offset, or A-offset, where A is the coordinate
        map(
            pair(parse_i32, opt(pair(alt((char('+'), char('-'))), parse_i32))),
            |(coordinate, offset)| {
                let offset = offset
                    .map(|(sign, value)| if sign == '-' { -value } else { value })
                    .unwrap_or(0);

                NucleotideCoordinate {
                    anchor: NucleotideAnchor::Absolute,
                    coordinate: CoordinateKind::Known(coordinate),
                    offset,
                }
            },
        ),
        // ?, unknown coordinate
        value(
            NucleotideCoordinate {
                anchor: NucleotideAnchor::Absolute,
                coordinate: CoordinateKind::Unknown,
                offset: 0,
            },
            char('?'),
        ),
    ))
    .parse(input)
}

/// Parses an unsigned decimal integer into `i32`.
fn parse_i32(input: &str) -> ParseResult<'_, i32> {
    map_res(digit1, str::parse::<i32>).parse(input)
}

fn parse_nonzero_i32(input: &str) -> ParseResult<'_, i32> {
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
        nucleotide_repeat_without_sequence,
        nucleotide_repeat_with_sequence,
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

/// Parses repeated-sequence syntax where only copy counts are written after
/// the top-level location.
fn nucleotide_repeat_without_sequence(input: &str) -> ParseResult<'_, NucleotideEdit> {
    map(
        pair(repeat_count_only_block, many0(repeat_located_count_block)),
        |(first, rest)| {
            let mut blocks = Vec::with_capacity(rest.len() + 1);
            blocks.push(first);
            blocks.extend(rest);
            NucleotideEdit::Repeat { blocks }
        },
    )
    .parse(input)
}

/// Parses repeated-sequence syntax where the repeated unit is written explicitly.
fn nucleotide_repeat_with_sequence(input: &str) -> ParseResult<'_, NucleotideEdit> {
    map(many1(repeat_sequence_block), |blocks| {
        NucleotideEdit::Repeat { blocks }
    })
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

/// Parses one top-level repeat block with an explicit sequence unit.
fn repeat_sequence_block(input: &str) -> ParseResult<'_, NucleotideRepeatBlock> {
    map(
        pair(
            nucleotide_literal,
            delimited(char('['), parse_usize, char(']')),
        ),
        |(unit, count)| NucleotideRepeatBlock {
            count,
            unit: Some(unit),
            location: None,
        },
    )
    .parse(input)
}

/// Parses the first count-only top-level repeat block after the main location.
fn repeat_count_only_block(input: &str) -> ParseResult<'_, NucleotideRepeatBlock> {
    map(delimited(char('['), parse_usize, char(']')), |count| {
        NucleotideRepeatBlock {
            count,
            unit: None,
            location: None,
        }
    })
    .parse(input)
}

/// Parses an additional located count-only repeat block used by composite repeats.
fn repeat_located_count_block(input: &str) -> ParseResult<'_, NucleotideRepeatBlock> {
    map(
        pair(
            nucleotide_interval,
            delimited(char('['), parse_usize, char(']')),
        ),
        |(location, count)| NucleotideRepeatBlock {
            count,
            unit: None,
            location: Some(location),
        },
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

/// Parser for protein variant and allele description.
fn protein_description(input: &str) -> ParseResult<'_, VariantDescription> {
    alt((
        // Protein allele description with known consequence: [Ser68Arg;Asn594del]
        // Composed of protein variant descriptions
        |input| {
            allele_description(
                input,
                protein_initial_allele,
                next_protein_allele,
                VariantDescription::ProteinAllele,
            )
        },
        // Predicted allele description: (Ser68Arg)
        map(
            delimited(char('('), protein_known_consequence, char(')')),
            |effect| {
                VariantDescription::Protein(ProteinVariant {
                    is_predicted: true,
                    effect,
                })
            },
        ),
        // Protein variant description: Ser68Arg, 0, ?
        map(protein_effect, |effect| {
            VariantDescription::Protein(ProteinVariant {
                is_predicted: false,
                effect,
            })
        }),
    ))
    .parse(input)
}

/// Parser for one protein variant description.
///
/// This parser admits three cases:
/// - No protein produced
/// - Known protein consequence, such as `Ser68Arg` or `0`
/// - Predicted protein consequence, such as `(Ser68Arg)`
///
/// One remaining gap: it does not admit plain `?`, which HGVS uses inside
/// protein alleles for an unknown allele member rather than a supported
/// protein consequence.
fn protein_variant_description(input: &str) -> ParseResult<'_, ProteinVariant> {
    alt((
        map(
            delimited(char('('), protein_known_consequence, char(')')),
            |effect| ProteinVariant {
                is_predicted: true,
                effect,
            },
        ),
        map(
            value(ProteinEffect::NoProteinProduced, char('0')),
            |effect| ProteinVariant {
                is_predicted: false,
                effect,
            },
        ),
        map(protein_known_consequence, |effect| ProteinVariant {
            is_predicted: false,
            effect,
        }),
    ))
    .parse(input)
}

/// Parser for digesting the initial allele written in a protein allele, such as:
///
/// - `p.[Ser68Arg;Asn594del]`
/// - `p.(Ser73Arg)(;)(Asn103del)`
fn protein_initial_allele(input: &str) -> ParseResult<'_, (Allele<ProteinVariant>, bool)> {
    alt((
        map(
            alt((
                // [(Ser73Arg;Asn103del)]
                delimited(char('['), predicted_protein_allele, char(']')),
                // [Ser68Arg;Asn594del]
                protein_bracket_allele,
            )),
            // true: initial allele wrapped in bracket
            |allele| (allele, true),
        ),
        // Ser68Arg
        map(protein_variant_description, |variant| {
            // false: initial allele not wrapped in bracket
            (Allele::from_variants(vec![variant]), false)
        }),
    ))
    .parse(input)
}

/// Parser for pattern of protein allele wrapped in bracket, such as:
///
/// - `[Ser68Arg;Asn594del]`
fn protein_bracket_allele(input: &str) -> ParseResult<'_, Allele<ProteinVariant>> {
    map(
        delimited(
            char('['),
            |input| variants_on_allele(input, protein_variant_description),
            char(']'),
        ),
        Allele::from_variants,
    )
    .parse(input)
}

/// Parser for the predicted protein allele pattern, such as:
///
/// - `(Ser73Arg;Asn103del)`.
fn predicted_protein_allele(input: &str) -> ParseResult<'_, Allele<ProteinVariant>> {
    map(
        delimited(
            char('('),
            |input| variants_on_allele(input, protein_known_consequence),
            char(')'),
        ),
        |effects| {
            Allele::from_variants(
                effects
                    .into_iter()
                    .map(|effect| ProteinVariant {
                        is_predicted: true,
                        effect,
                    })
                    .collect(),
            )
        },
    )
    .parse(input)
}

/// Parser for digesting the next protein allele following the initial allele +
/// phase marker in the protein allele description, such as:
///
/// - `[Ser68_Arg70=]` in `p.[Ser68_Arg70dup];[Ser68_Arg70=]` after `;`
///
/// Like in the nucleotide case, the phase marker after the initial allele
/// imposes syntactic rule.
fn next_protein_allele(phase: AllelePhase, input: &str) -> ParseResult<'_, Allele<ProteinVariant>> {
    match phase {
        AllelePhase::Trans => alt((
            // [(Ser73Arg;Asn103del)]
            delimited(char('['), predicted_protein_allele, char(']')),
            // [Ser68Arg;Asn594del]
            protein_bracket_allele,
        ))
        .parse(input),
        AllelePhase::Uncertain => map(protein_variant_description, |variant| {
            Allele::from_variants(vec![variant])
        })
        .parse(input),
    }
}

/// Parser for one known protein consequence, which means:
///
/// - `Ser68Arg`
/// - `Ser68del`
/// - `Ser68_Ala74insSerGln`
///
/// The parser does not digest `?` and `0`.
fn protein_known_consequence(input: &str) -> ParseResult<'_, ProteinEffect> {
    map_res(
        pair(protein_location, protein_edit),
        build_protein_edit_effect,
    )
    .parse(input)
}

/// Parses the supported protein effect types. Differentiating unknown consequence
/// from the case where mutation produces no protein.
fn protein_effect(input: &str) -> ParseResult<'_, ProteinEffect> {
    alt((
        value(ProteinEffect::Unknown, char('?')),
        value(ProteinEffect::NoProteinProduced, char('0')),
        protein_known_consequence,
    ))
    .parse(input)
}

fn build_protein_edit_effect(
    (location, edit): (Location<ProteinCoordinate>, ProteinEdit),
) -> Result<ProteinEffect, ()> {
    let location = resolve_protein_effect_location(&location, &edit).ok_or(())?;
    Ok(ProteinEffect::Edit { location, edit })
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
        map(
            separated_pair(protein_coordinate, char('_'), protein_coordinate),
            |(start, end)| Interval {
                start,
                end: Some(end),
            },
        ),
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
        map(
            separated_pair(
                protein_uncertain_interval,
                char('_'),
                protein_uncertain_interval,
            ),
            |(start, end)| Interval {
                start,
                end: Some(end),
            },
        ),
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
        map(delimited(char('['), parse_usize, char(']')), |count| {
            ProteinEdit::Repeat { count }
        }),
        protein_extension_edit,
        protein_frameshift_edit,
        map(preceded(tag("ins"), protein_sequence), |sequence| {
            ProteinEdit::Insertion { sequence }
        }),
        map(protein_symbol, |to| ProteinEdit::Substitution { to }),
    ))
    .parse(input)
}

/// Parses N-terminal and C-terminal protein extension syntax.
fn protein_extension_edit(input: &str) -> ParseResult<'_, ProteinEdit> {
    alt((
        map(
            preceded(tag("ext"), protein_n_terminal_extension_ordinal),
            |terminal_ordinal| {
                ProteinEdit::Extension(ProteinExtensionEdit {
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
                ProteinEdit::Extension(ProteinExtensionEdit {
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
fn protein_frameshift_edit(input: &str) -> ParseResult<'_, ProteinEdit> {
    alt((
        map(
            pair(
                protein_frameshift_residue,
                pair(tag("fs"), protein_frameshift_stop),
            ),
            |(to_residue, (_, stop))| ProteinEdit::Frameshift {
                to_residue: Some(to_residue),
                stop,
            },
        ),
        value(
            ProteinEdit::Frameshift {
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
            preceded(alt((tag("Ter"), tag("*"))), parse_usize),
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

fn resolve_protein_effect_location(
    location: &Location<ProteinCoordinate>,
    edit: &ProteinEdit,
) -> Option<Location<ProteinCoordinate>> {
    let ProteinEdit::Extension(extension) = edit else {
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

/// Returns the top-level nucleotide location to expose on the parsed variant.
fn resolve_nucleotide_location(
    initial_location: &Location<NucleotideCoordinate>,
    edit: &NucleotideEdit,
) -> Option<Location<NucleotideCoordinate>> {
    let NucleotideEdit::Repeat { blocks } = edit else {
        return Some(initial_location.clone());
    };

    let Location::Known(initial_location) = initial_location else {
        return None;
    };

    let Some(last_location) = blocks
        .iter()
        .filter_map(|block| block.location.as_ref())
        .last()
    else {
        return Some(Location::from_known(initial_location.clone()));
    };

    Some(Location::from_known(Interval {
        start: initial_location.start.clone(),
        end: last_location
            .end
            .clone()
            .or_else(|| Some(last_location.start.clone())),
    }))
}

/// Validates molecule-specific rules for repeated-sequence descriptions.
fn is_valid_nucleotide_repeat(
    coordinate_system: CoordinateSystem,
    initial_location: &Location<NucleotideCoordinate>,
    edit: &NucleotideEdit,
) -> bool {
    let NucleotideEdit::Repeat { blocks } = edit else {
        return true;
    };

    let Location::Known(initial_location) = initial_location else {
        return false;
    };

    let all_have_units = blocks.iter().all(|block| block.unit.is_some());
    let none_have_units = blocks.iter().all(|block| block.unit.is_none());
    let any_have_locations = blocks.iter().any(|block| block.location.is_some());

    match coordinate_system {
        CoordinateSystem::Rna => {
            if none_have_units {
                true
            } else if all_have_units {
                blocks.len() == 1 && !any_have_locations && initial_location.end.is_none()
            } else {
                false
            }
        }
        CoordinateSystem::Genomic
        | CoordinateSystem::CircularGenomic
        | CoordinateSystem::Mitochondrial
        | CoordinateSystem::CodingDna
        | CoordinateSystem::NonCodingDna => all_have_units && !any_have_locations,
        CoordinateSystem::Protein => false,
    }
}

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;

    use super::*;

    #[test]
    fn parses_nucleotide_position_branches() {
        let (_, coding) = all_consuming(nucleotide_coordinate).parse("93+1").unwrap();
        assert_eq!(coding.anchor, NucleotideAnchor::Absolute);
        assert_eq!(coding.coordinate, CoordinateKind::Known(93));
        assert_eq!(coding.offset, 1);

        let (_, upstream_intronic) = all_consuming(nucleotide_coordinate).parse("93-2").unwrap();
        assert_eq!(upstream_intronic.anchor, NucleotideAnchor::Absolute);
        assert_eq!(upstream_intronic.coordinate, CoordinateKind::Known(93));
        assert_eq!(upstream_intronic.offset, -2);

        let (_, utr5) = all_consuming(nucleotide_coordinate).parse("-18").unwrap();
        assert_eq!(utr5.anchor, NucleotideAnchor::RelativeCdsStart);
        assert_eq!(utr5.coordinate, CoordinateKind::Known(-18));
        assert_eq!(utr5.offset, 0);

        let (_, utr5_intronic) = all_consuming(nucleotide_coordinate)
            .parse("-106+2")
            .unwrap();
        assert_eq!(utr5_intronic.anchor, NucleotideAnchor::RelativeCdsStart);
        assert_eq!(utr5_intronic.coordinate, CoordinateKind::Known(-106));
        assert_eq!(utr5_intronic.offset, 2);

        let (_, utr5_intronic_upstream) =
            all_consuming(nucleotide_coordinate).parse("-84-1").unwrap();
        assert_eq!(
            utr5_intronic_upstream.anchor,
            NucleotideAnchor::RelativeCdsStart
        );
        assert_eq!(
            utr5_intronic_upstream.coordinate,
            CoordinateKind::Known(-84)
        );
        assert_eq!(utr5_intronic_upstream.offset, -1);

        let (_, utr3) = all_consuming(nucleotide_coordinate).parse("*18").unwrap();
        assert_eq!(utr3.anchor, NucleotideAnchor::RelativeCdsEnd);
        assert_eq!(utr3.coordinate, CoordinateKind::Known(18));
        assert_eq!(utr3.offset, 0);

        let (_, utr3_intronic) = all_consuming(nucleotide_coordinate)
            .parse("*639-1")
            .unwrap();
        assert_eq!(utr3_intronic.anchor, NucleotideAnchor::RelativeCdsEnd);
        assert_eq!(utr3_intronic.coordinate, CoordinateKind::Known(639));
        assert_eq!(utr3_intronic.offset, -1);

        assert!(all_consuming(nucleotide_coordinate).parse("-0").is_err());
        assert!(all_consuming(nucleotide_coordinate).parse("-0+2").is_err());
        assert!(all_consuming(nucleotide_coordinate).parse("*0").is_err());
        assert!(all_consuming(nucleotide_coordinate).parse("*0-1").is_err());
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
        assert!(matches!(
            all_consuming(nucleotide_edit).parse("[4]").unwrap().1,
            NucleotideEdit::Repeat { .. }
        ));
        assert!(matches!(
            all_consuming(nucleotide_edit).parse("CAG[23]").unwrap().1,
            NucleotideEdit::Repeat { .. }
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
        assert!(matches!(
            all_consuming(protein_effect).parse("Ala2[10]").unwrap().1,
            ProteinEffect::Edit {
                edit: ProteinEdit::Repeat { count: 10 },
                ..
            }
        ));
        assert!(matches!(
            all_consuming(protein_effect).parse("Arg97fs").unwrap().1,
            ProteinEffect::Edit {
                edit: ProteinEdit::Frameshift {
                    to_residue: None,
                    stop: ProteinFrameshiftStop {
                        ordinal: None,
                        kind: ProteinFrameshiftStopKind::Omitted,
                    },
                },
                ..
            }
        ));
        assert!(matches!(
            all_consuming(protein_effect).parse("Met1ext-5").unwrap().1,
            ProteinEffect::Edit {
                edit: ProteinEdit::Extension(ProteinExtensionEdit {
                    to_terminal: ProteinExtensionTerminal::N,
                    to_residue: None,
                    terminal_ordinal: Some(-5),
                }),
                ..
            }
        ));
        assert!(matches!(
            all_consuming(protein_effect)
                .parse("Ter110GlnextTer17")
                .unwrap()
                .1,
            ProteinEffect::Edit {
                edit: ProteinEdit::Extension(ProteinExtensionEdit {
                    to_terminal: ProteinExtensionTerminal::C,
                    to_residue: Some(_),
                    terminal_ordinal: Some(17),
                }),
                ..
            }
        ));
        assert!(matches!(
            all_consuming(protein_effect)
                .parse("Arg97ProfsTer23")
                .unwrap()
                .1,
            ProteinEffect::Edit {
                edit: ProteinEdit::Frameshift {
                    to_residue: Some(_),
                    stop: ProteinFrameshiftStop {
                        ordinal: Some(23),
                        kind: ProteinFrameshiftStopKind::Known,
                    },
                },
                ..
            }
        ));
        assert!(matches!(
            all_consuming(protein_effect)
                .parse("Ile327Argfs*?")
                .unwrap()
                .1,
            ProteinEffect::Edit {
                edit: ProteinEdit::Frameshift {
                    to_residue: Some(_),
                    stop: ProteinFrameshiftStop {
                        ordinal: None,
                        kind: ProteinFrameshiftStopKind::Unknown,
                    },
                },
                ..
            }
        ));
    }

    #[test]
    fn parses_protein_frameshift_branches() {
        assert_eq!(
            all_consuming(protein_edit).parse("fs").unwrap().1,
            ProteinEdit::Frameshift {
                to_residue: None,
                stop: ProteinFrameshiftStop {
                    ordinal: None,
                    kind: ProteinFrameshiftStopKind::Omitted,
                },
            }
        );
        assert_eq!(
            all_consuming(protein_edit).parse("ProfsTer23").unwrap().1,
            ProteinEdit::Frameshift {
                to_residue: Some("Pro".to_string()),
                stop: ProteinFrameshiftStop {
                    ordinal: Some(23),
                    kind: ProteinFrameshiftStopKind::Known,
                },
            }
        );
        assert_eq!(
            all_consuming(protein_edit).parse("Argfs*?").unwrap().1,
            ProteinEdit::Frameshift {
                to_residue: Some("Arg".to_string()),
                stop: ProteinFrameshiftStop {
                    ordinal: None,
                    kind: ProteinFrameshiftStopKind::Unknown,
                },
            }
        );
        assert!(all_consuming(protein_edit).parse("TerfsTer2").is_err());
    }

    #[test]
    fn parses_protein_extension_branches() {
        assert_eq!(
            all_consuming(protein_edit).parse("ext-5").unwrap().1,
            ProteinEdit::Extension(ProteinExtensionEdit {
                to_terminal: ProteinExtensionTerminal::N,
                to_residue: None,
                terminal_ordinal: Some(-5),
            })
        );
        assert_eq!(
            all_consuming(protein_edit).parse("GlnextTer17").unwrap().1,
            ProteinEdit::Extension(ProteinExtensionEdit {
                to_terminal: ProteinExtensionTerminal::C,
                to_residue: Some("Gln".to_string()),
                terminal_ordinal: Some(17),
            })
        );
        assert_eq!(
            all_consuming(protein_edit).parse("Argext*?").unwrap().1,
            ProteinEdit::Extension(ProteinExtensionEdit {
                to_terminal: ProteinExtensionTerminal::C,
                to_residue: Some("Arg".to_string()),
                terminal_ordinal: None,
            })
        );
        assert!(all_consuming(protein_edit).parse("TerextTer17").is_err());
    }
}
