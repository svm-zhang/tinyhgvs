//! Lightweight post-parse diagnostics for unsupported HGVS families.
//!
//! The parser intentionally recognizes only supported syntax. When parsing
//! fails, this module runs a second, shallow classification pass over the raw
//! input to produce stable diagnostic codes such as `unsupported.telomeric_position`.
//!
//! This is not a second full parser. The matchers stay intentionally small and
//! ordered so they can be retired as syntax becomes truly supported.

use crate::error::ParseHgvsError;

/// One ordered matcher in the unsupported-syntax classifier.
struct DiagnosticMatcher {
    code: &'static str,
    message: &'static str,
    detect: fn(&str) -> Option<String>,
}

// Certain unsupported families should be captured before broader
// fallback cases.
const UNSUPPORTED_MATCHERS: &[DiagnosticMatcher] = &[
    // Example: `NM_002354.2:r.-358_555::NM_000251.2:r.212_*279`
    DiagnosticMatcher {
        code: "unsupported.rna_adjoined_transcript",
        message: "RNA adjoined transcript syntax is not supported yet",
        detect: rna_adjoined_transcript_fragment,
    },
    // Example: `NC_000023.11(NM_004006.2):r.[897u>g,832_960del]`
    DiagnosticMatcher {
        code: "unsupported.rna_splicing_outcome",
        message: "RNA splicing outcome containers are not supported yet",
        detect: rna_splicing_outcome_fragment,
    },
    // Example: `NC_000023.11:g.pter_qtersup`
    DiagnosticMatcher {
        code: "unsupported.telomeric_position",
        message: "telomeric positions such as pter and qter are not supported yet",
        detect: telomeric_position_fragment,
    },
    // Example: `NC_000011.10:g.1999904_1999946|gom`
    DiagnosticMatcher {
        code: "unsupported.epigenetic_edit",
        message: "epigenetic edit syntax is not supported yet",
        detect: epigenetic_edit_fragment,
    },
    // Examples: `p.Arg78_Gly79insXaa[23]`, `...ins*63`
    DiagnosticMatcher {
        code: "unsupported.protein_insertion_content",
        message: "quantified or terminal protein insertion content is not supported yet",
        detect: protein_insertion_content_fragment,
    },
    // Example: `p.(Gly719Ala^Ser)`
    DiagnosticMatcher {
        code: "unsupported.protein_uncertain_consequence",
        message: "uncertain protein consequence syntax is not supported yet",
        detect: protein_uncertain_consequence_fragment,
    },
];

/// Classifies a parse failure into a stable diagnostic code when possible.
pub(crate) fn classify_parse_failure(input: &str) -> ParseHgvsError {
    for matcher in UNSUPPORTED_MATCHERS {
        if let Some(fragment) = (matcher.detect)(input) {
            return ParseHgvsError::unsupported(
                matcher.code,
                matcher.message,
                input,
                Some(fragment),
            );
        }
    }

    ParseHgvsError::invalid(input)
}

/// Detects RNA adjoined transcript syntax such as `r.-358_555::r.212_*279`.
fn rna_adjoined_transcript_fragment(input: &str) -> Option<String> {
    let description = coordinate_description_fragment(input, "r.")?;
    (description.contains("::") || input.contains("::")).then(|| "::".to_string())
}

// These are higher-level transcript consequence containers. Direct RNA states
// such as `r.?`, `r.0`, and `r.spl` are parsed by the core parser now.
/// Detects still-unsupported top-level RNA splicing outcome containers.
fn rna_splicing_outcome_fragment(input: &str) -> Option<String> {
    let description = coordinate_description_fragment(input, "r.")?;
    let has_context = input.contains("):r.");

    if !has_context {
        return None;
    }

    if description.starts_with('[') && description.contains(',') {
        Some("r.[...]".to_string())
    } else if description.starts_with('(') && description.contains(',') {
        Some("r.(...)".to_string())
    } else {
        None
    }
}

/// Detects telomeric positions such as `pter` and `qter`.
fn telomeric_position_fragment(input: &str) -> Option<String> {
    if input.contains("pter") {
        Some("pter".to_string())
    } else if input.contains("qter") {
        Some("qter".to_string())
    } else {
        None
    }
}

/// Detects epigenetic edit modifiers introduced with `|`.
fn epigenetic_edit_fragment(input: &str) -> Option<String> {
    let description = variant_description_fragment(input)?;
    description
        .split_once('|')
        .map(|(_, modifier)| format!("|{modifier}"))
}

/// Detects unsupported quantified or terminal protein insertion content.
fn protein_insertion_content_fragment(input: &str) -> Option<String> {
    let description = protein_description_fragment(input)?;
    if !description.contains("ins") {
        return None;
    }

    if description.contains("Xaa[") {
        Some("Xaa[...]".to_string())
    } else if description.contains("ins*") {
        Some("*".to_string())
    } else {
        None
    }
}

/// Detects uncertain protein consequences written with `^`.
fn protein_uncertain_consequence_fragment(input: &str) -> Option<String> {
    let description = protein_description_fragment(input)?;
    if description.starts_with('[') {
        return None;
    }

    if description.contains('^') {
        Some("^".to_string())
    } else if description.contains("[(") {
        Some("[(...)]".to_string())
    } else {
        None
    }
}

fn protein_description_fragment(input: &str) -> Option<&str> {
    coordinate_description_fragment(input, "p.")
}

fn variant_description_fragment(input: &str) -> Option<&str> {
    ["g.", "o.", "m.", "c.", "n.", "r.", "p."]
        .into_iter()
        .find_map(|marker| coordinate_description_fragment(input, marker))
}

fn coordinate_description_fragment<'a>(input: &'a str, marker: &str) -> Option<&'a str> {
    if let Some(rest) = input.strip_prefix(marker) {
        return Some(rest);
    }

    let needle = format!(":{marker}");
    input
        .find(&needle)
        .map(|index| &input[index + needle.len()..])
}
