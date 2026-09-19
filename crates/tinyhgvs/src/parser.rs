//! HGVS variant parsers.
//!
//! The parsing strategy:
//!
//! - parse the accepted syntax into the Rust data model
//! - route rejected inputs to the lightweight diagnostic classifier

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::char;
use nom::combinator::{all_consuming, map, opt};
use nom::sequence::{pair, terminated};
use nom::{IResult, Parser};

use crate::diagnostics::classify_parse_failure;
use crate::error::ParseHgvsError;
use crate::model::{CoordinateSystem, HgvsVariant};

type ParseResult<'a, T> = IResult<&'a str, T>;

mod cdna;
mod core;
mod location;
mod nucleotide;
mod protein;
mod repeat;
mod rna;

use cdna::cdna_description;
use core::reference_spec;
use nucleotide::genomic_description;
use protein::protein_description;
use rna::rna_description;

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
/// let error = parse_hgvs("NC_000023.11:g.pter_qtersup").unwrap_err();
/// assert_eq!(error.code(), "unsupported.telomeric_position");
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
