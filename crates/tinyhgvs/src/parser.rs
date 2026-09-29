//! Public parser entry point and private grammar modules.
//!
//! [`parse_hgvs`] is the crate-facing parser. The submodules below organize
//! private grammar families for references, locations, repeats, nucleotide
//! descriptions, RNA descriptions, and protein descriptions.

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
/// A coding-DNA substitution with an intronic offset:
///
/// ```rust
/// use tinyhgvs::{
///     CodingDnaOutcome, NucleotideAnchor, NucleotideEditKind, VariantDescription, parse_hgvs,
/// };
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("  NM_004006.2:c.357+1G>A  ")?;
///
/// let VariantDescription::CodingDna(CodingDnaOutcome::Known(edit)) = variant.description else {
///     panic!("expected a coding-DNA edit");
/// };
///
/// assert_eq!(edit.location.start().unwrap().anchor(), Some(NucleotideAnchor::Absolute));
/// assert_eq!(edit.location.start().unwrap().coordinate(), Some(357));
/// assert_eq!(edit.location.start().unwrap().offset(), Some(1));
/// assert!(matches!(
///     edit.kind,
///     NucleotideEditKind::Substitution { ref reference, ref alternate }
///         if reference == "G" && alternate == "A"
/// ));
/// # Ok(())
/// # }
/// ```
///
/// A coding-DNA allele variant with two in-trans alleles:
///
/// ```rust
/// use tinyhgvs::{AlleleForm, AllelePhase, VariantDescription, parse_hgvs};
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("NM_004006.2:c.[2376G>C];[2376=]")?;
///
/// let VariantDescription::CodingDnaAllele(AlleleForm::Single(allele)) = variant.description else {
///     panic!("expected a coding-DNA allele");
/// };
///
/// assert_eq!(allele.phase, Some(AllelePhase::Trans));
/// assert_eq!(allele.allele_one.variants.len(), 1);
/// assert_eq!(allele.allele_two.as_ref().unwrap().variants.len(), 1);
/// # Ok(())
/// # }
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
