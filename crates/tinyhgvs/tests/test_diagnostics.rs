use tinyhgvs::{parse_hgvs, ParseHgvsErrorKind};

fn parse_error(example: &str) -> tinyhgvs::ParseHgvsError {
    parse_hgvs(example).unwrap_err()
}

#[test]
fn classifies_supported_diagnostic_codes() {
    let cases = [
        (
            "NC_000023.11:g.pter_qtersup",
            "unsupported.telomeric_position",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "telomeric positions such as pter and qter are not supported yet",
            Some("pter"),
        ),
        (
            "NC_000011.10:g.1999904_1999946|gom",
            "unsupported.epigenetic_edit",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "epigenetic edit syntax is not supported yet",
            Some("|gom"),
        ),
        (
            "NM_004006.3:r.spl",
            "unsupported.rna_special_state",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "RNA consequence states such as r.spl, r.?, and r.0 are not supported yet",
            Some("r.spl"),
        ),
        (
            "NC_000023.11(NM_004006.2):r.[897u>g,832_960del]",
            "unsupported.rna_splicing_outcome",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "RNA splicing outcome containers are not supported yet",
            Some("r.[...]"),
        ),
        (
            "NM_002354.2:r.-358_555::NM_000251.2:r.212_*279",
            "unsupported.rna_adjoined_transcript",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "RNA adjoined transcript syntax is not supported yet",
            Some("::"),
        ),
        (
            "NM_004006.2:c.[2376G>C];[?]",
            "unsupported.allele_unknown_variant",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "allele variants written as [?] are not supported yet",
            Some("[?]"),
        ),
        (
            "NM_004006.2:c.2376G>C(;)(2376G>C)",
            "unsupported.allele_uncertain_variant_state",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "uncertain allele variant states are not supported yet",
            Some("(;)(...)"),
        ),
        (
            "r.-124_-123[14];[18]",
            "unsupported.allele",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "allele syntax is not supported yet",
            Some("];["),
        ),
        (
            "NP_003997.1:p.[Lys31Asn,Val25_Lys31del]",
            "unsupported.one_allele_multi_protein",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "one protein allele encoding more than one protein is not supported yet",
            Some(","),
        ),
        (
            "NP_003997.2:p.[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]",
            "unsupported.alternate_allele_state",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "alternate allele states are not supported yet",
            Some("^"),
        ),
        (
            "r.-128_-126[(600_800)]",
            "unsupported.uncertain_size",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "uncertain HGVS size syntax is not supported yet",
            Some("[(...)]"),
        ),
        (
            "p.Arg78_Gly79insXaa[23]",
            "unsupported.protein_insertion_payload",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "quantified or terminal protein insertion payloads are not supported yet",
            Some("Xaa[...]"),
        ),
        (
            "p.(Gly719Ala^Ser)",
            "unsupported.protein_uncertain_consequence",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "uncertain protein consequence syntax is not supported yet",
            Some("^"),
        ),
        (
            "p.(Gln18)[(70_80)]",
            "unsupported.protein_uncertain_consequence",
            ParseHgvsErrorKind::UnsupportedSyntax,
            "uncertain protein consequence syntax is not supported yet",
            Some("[(...)]"),
        ),
    ];

    for (input, code, kind, message, fragment) in cases {
        let error = parse_error(input);
        assert_eq!(error.code(), code, "unexpected code for {input}");
        assert_eq!(error.kind(), kind, "unexpected kind for {input}");
        assert_eq!(error.message(), message, "unexpected message for {input}");
        assert_eq!(error.input(), input, "unexpected input for {input}");
        assert_eq!(
            error.fragment(),
            fragment,
            "unexpected fragment for {input}"
        );
        assert_eq!(
            error.parser_version(),
            env!("CARGO_PKG_VERSION"),
            "unexpected parser version for {input}"
        );
    }
}

#[test]
fn prioritizes_specific_rna_codes_before_generic_ones() {
    let splicing = parse_error("NC_000023.11(NM_004006.2):r.spl");
    let unknown_member = parse_error("NM_004006.2:c.[2376G>C];[?]");
    let protein_unknown_member = parse_error("NP_003997.1:p.[(Ser68Arg)];[?]");
    let uncertain_state = parse_error("NM_004006.2:c.2376G>C(;)(2376G>C)");

    assert_eq!(splicing.code(), "unsupported.rna_splicing_outcome");
    assert_eq!(unknown_member.code(), "unsupported.allele_unknown_variant");
    assert_eq!(
        protein_unknown_member.code(),
        "unsupported.allele_unknown_variant"
    );
    assert_eq!(
        uncertain_state.code(),
        "unsupported.allele_uncertain_variant_state"
    );
}

#[test]
fn falls_back_to_invalid_syntax_for_unclassified_failures() {
    let error = parse_error("not-hgvs");

    assert_eq!(error.code(), "invalid.syntax");
    assert_eq!(error.kind(), ParseHgvsErrorKind::InvalidSyntax);
    assert_eq!(error.message(), "failed to parse HGVS variant");
    assert_eq!(error.input(), "not-hgvs");
    assert_eq!(error.fragment(), None);
}

#[test]
fn malformed_uncertain_range_now_falls_back_to_generic_invalid_syntax() {
    let cases = [
        "NC_000023.10:g.(?_?)del",
        "NC_000023.10:g.(?_?)_(?_?)del",
        "NM_004006.2:r.(?_?)del",
        "NM_004006.2:r.(?_?)_(?_?)del",
    ];

    for input in cases {
        let error = parse_error(input);

        assert_eq!(error.code(), "invalid.syntax");
        assert_eq!(error.kind(), ParseHgvsErrorKind::InvalidSyntax);
        assert_eq!(error.message(), "failed to parse HGVS variant");
        assert_eq!(error.fragment(), None);
    }
}

#[test]
fn displays_machine_code_message_and_version() {
    let error = parse_error("p.Arg78_Gly79insXaa[23]");
    let rendered = error.to_string();

    assert!(rendered.contains("[unsupported.protein_insertion_payload]"));
    assert!(rendered
        .contains("quantified or terminal protein insertion payloads are not supported yet"));
    assert!(rendered.contains("`p.Arg78_Gly79insXaa[23]`"));
    assert!(rendered.contains(env!("CARGO_PKG_VERSION")));
}
