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
    // Examples: `NM_001385026.1:c.-666+629C>T`, `...:c.*24-12888C>T`
    DiagnosticMatcher {
        code: "unsupported.cdna_offset_anchor",
        message: "coding-DNA positions anchored to CDS start/end with additional offsets are not supported yet",
        detect: cdna_offset_anchor_fragment,
    },
    // Examples: `g.(?_234567)_(345678_?)del`, `g.(33038277_33038278)C>T`
    DiagnosticMatcher {
        code: "unsupported.uncertain_range",
        message: "uncertain HGVS ranges are not supported yet",
        detect: uncertain_range_fragment,
    },
    // Examples: `g.[123G>A;345del]`, `r.[123c>a;345del]`, `p.Val7=/del`
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

// This catches CDS-start/CDS-end anchored positions with additional offsets,
// which are valid HGVS patterns but not representable by the current position type.
/// Detects unsupported CDS-start or CDS-end anchored offsets.
fn cdna_offset_anchor_fragment(input: &str) -> Option<String> {
    let description = coordinate_description_fragment(input, "c.")?;
    let first = description
        .split_once('_')
        .map_or(description, |(start, _)| start);
    if let Some(fragment) = anchored_offset_position_fragment(first) {
        return Some(fragment);
    }

    description
        .split_once('_')
        .and_then(|(_, end)| anchored_offset_position_fragment(end))
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

// A single top-level code is used across DNA, RNA, and protein until the model
// grows explicit allele and phase containers.
/// Detects allele containers across DNA, RNA, and protein syntax.
fn allele_fragment(input: &str) -> Option<String> {
    let description = variant_description_fragment(input)?;

    if description.contains("=/") {
        Some("=/".to_string())
    } else if description.contains("(;)") {
        Some("(;)".to_string())
    } else if description.contains("];[") {
        Some("];[".to_string())
    } else if description.starts_with('[')
        && (description.contains(';') || description.contains(','))
    {
        Some("[".to_string())
    } else {
        None
    }
}

/// Extracts the offset fragment from a CDS-anchored position like `-666+629`.
fn anchored_offset_position_fragment(position: &str) -> Option<String> {
    let mut characters = position.char_indices();
    let (_, anchor) = characters.next()?;

    if !matches!(anchor, '-' | '*') {
        return None;
    }

    let mut offset_start = None;

    for (index, character) in characters.by_ref() {
        if character.is_ascii_digit() {
            continue;
        }

        if matches!(character, '+' | '-') {
            offset_start = Some(index);
        }
        break;
    }

    let offset_start = offset_start?;

    let mut end = offset_start + 1;
    for (index, character) in characters {
        if character.is_ascii_digit() {
            end = index + character.len_utf8();
            continue;
        }
        break;
    }

    (end > offset_start + 1).then(|| position[..end].to_string())
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
