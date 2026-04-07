//! Lightweight post-parse diagnostics for unsupported HGVS families.
//!
//! The parser intentionally recognizes only supported syntax. When parsing
//! fails, this module runs a second, shallow classification pass over the raw
//! input to produce stable diagnostic codes such as `unsupported.allele`.
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
    // Examples: `NC_000023.11(NM_004006.2):r.[...]`, `...:r.spl`
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
        code: "unsupported.protein_insertion_payload",
        message: "quantified or terminal protein insertion payloads are not supported yet",
        detect: protein_insertion_payload_fragment,
    },
    // Example: `p.(Gly719Ala^Ser)`
    DiagnosticMatcher {
        code: "unsupported.protein_uncertain_consequence",
        message: "uncertain protein consequence syntax is not supported yet",
        detect: protein_uncertain_consequence_fragment,
    },
    // Example: `NM_004006.2:r.(222_226)insg`
    DiagnosticMatcher {
        code: "unsupported.rna_uncertain_position",
        message: "RNA variants with uncertain positions are not supported yet",
        detect: rna_uncertain_position_fragment,
    },
    // Examples: `NM_004006.3:r.spl`, `r.?`, `r.0`, `r.(1388g>a)`
    DiagnosticMatcher {
        code: "unsupported.rna_special_state",
        message: "RNA consequence states such as r.spl, r.?, and r.0 are not supported yet",
        detect: rna_special_state_fragment,
    },
    // Examples: `g.(?_234567)_(345678_?)del`, `g.(33038277_33038278)C>T`
    DiagnosticMatcher {
        code: "unsupported.uncertain_range",
        message: "uncertain HGVS ranges are not supported yet",
        detect: uncertain_range_fragment,
    },
    // Example: `c.[2376G>C];[?]`
    DiagnosticMatcher {
        code: "unsupported.allele_unknown_variant",
        message: "allele variants written as [?] are not supported yet",
        detect: allele_unknown_variant_fragment,
    },
    // Example: `c.2376G>C(;)(2376G>C)`
    DiagnosticMatcher {
        code: "unsupported.allele_uncertain_variant_state",
        message: "uncertain allele variant states are not supported yet",
        detect: allele_uncertain_variant_state_fragment,
    },
    // Examples: `p.Val7=/del`, `p.[(Ser73Arg;Asn103del)]`
    DiagnosticMatcher {
        code: "unsupported.protein_allele",
        message: "protein allele syntax is not supported yet",
        detect: protein_allele_fragment,
    },
    // Example: `r.-124_-123[14];[18]`
    DiagnosticMatcher {
        code: "unsupported.allele",
        message: "allele syntax is not supported yet",
        detect: allele_fragment,
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

// These are higher-level transcript consequence containers, not the exact
// nucleotide edits that the parser already supports for some splice outcomes.
/// Detects top-level RNA splicing outcome containers.
fn rna_splicing_outcome_fragment(input: &str) -> Option<String> {
    let description = coordinate_description_fragment(input, "r.")?;
    let has_context = input.contains("):r.");

    if !has_context {
        return None;
    }

    if description.starts_with('[') {
        Some("r.[...]".to_string())
    } else if description.starts_with('(') {
        Some("r.(...)".to_string())
    } else if description.starts_with('?') {
        Some("r.?".to_string())
    } else if description.starts_with("spl") {
        Some("r.spl".to_string())
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

/// Detects unsupported quantified or terminal protein insertion payloads.
fn protein_insertion_payload_fragment(input: &str) -> Option<String> {
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
    if description.contains('^') {
        Some("^".to_string())
    } else if description.contains("[(") {
        Some("[(...)]".to_string())
    } else {
        None
    }
}

/// Detects RNA edits with uncertain positions such as `r.(222_226)insg`.
fn rna_uncertain_position_fragment(input: &str) -> Option<String> {
    let description = coordinate_description_fragment(input, "r.")?;
    (description.starts_with('(')
        && (description.contains("ins")
            || description.contains("del")
            || description.contains("dup")
            || description.contains("inv")))
    .then(|| "r.(...)".to_string())
}

// The current RNA model assumes exact `location + edit`, so special states and
// slash-style outcome forms are grouped here until RNA consequence types expand.
/// Detects RNA special states such as `r.?`, `r.spl`, and `r.0`.
fn rna_special_state_fragment(input: &str) -> Option<String> {
    let description = coordinate_description_fragment(input, "r.")?;

    if description.starts_with('?') {
        Some("r.?".to_string())
    } else if description.starts_with("spl") {
        Some("r.spl".to_string())
    } else if description.starts_with('0') {
        Some("r.0".to_string())
    } else if description.starts_with('(') {
        Some("r.(...)".to_string())
    } else if description.contains("=/") || description.contains("//") {
        Some("=/".to_string())
    } else {
        None
    }
}

// This remains deliberately broad until uncertainty gets its own richer model.
/// Detects uncertainty wrappers and range forms handled by the broad fallback.
fn uncertain_range_fragment(input: &str) -> Option<String> {
    let description = variant_description_fragment(input)?;

    if description.contains('^') {
        Some("^".to_string())
    } else if description.contains("[(") {
        Some("[(...)]".to_string())
    } else if description.contains("(?") || description.contains("?)") {
        Some("(?)".to_string())
    } else if description.starts_with('(') || description.contains("_(") {
        Some("(".to_string())
    } else {
        None
    }
}

/// Detects allele variants written as `[?]`.
fn allele_unknown_variant_fragment(input: &str) -> Option<String> {
    let description = variant_description_fragment(input)?;
    description.contains("[?]").then(|| "[?]".to_string())
}

/// Detects allele forms where a variant is written but its allele state is uncertain.
fn allele_uncertain_variant_state_fragment(input: &str) -> Option<String> {
    let description = variant_description_fragment(input)?;
    description.contains("(;)(").then(|| "(;)(...)".to_string())
}

/// Detects protein allele containers, which remain out of scope.
fn protein_allele_fragment(input: &str) -> Option<String> {
    let description = protein_description_fragment(input)?;
    allele_like_fragment(description)
}

/// Detects other still-unsupported allele containers.
fn allele_fragment(input: &str) -> Option<String> {
    let description = variant_description_fragment(input)?;
    allele_like_fragment(description)
}

fn allele_like_fragment(description: &str) -> Option<String> {
    if description.contains("=/") {
        Some("=/".to_string())
    } else if description.contains("];[") {
        Some("];[".to_string())
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
