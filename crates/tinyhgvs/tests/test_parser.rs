mod utils;

use tinyhgvs::{
    parse_hgvs, AlleleForm, AllelePhase, AlleleStateCertainty, CoordinateSystem, NucleotideAnchor,
    NucleotideEditKind, NucleotideSequenceItem,
    OutcomeCertainty, ProteinEditKind, ProteinExtensionTerminal, ProteinFrameshiftStopKind,
    ProteinOutcome, RnaOutcome,
};
use utils::prelude::*;

#[test]
fn parses_reference_context() {
    let variant = parse_variant("NG_012232.1(NM_004006.2):c.93+1G>T");

    assert_eq!(variant.coordinate_system, CoordinateSystem::CodingDna);
    let reference = variant.reference.expect("expected reference");
    assert_eq!(reference.primary.id, "NG_012232.1");
    assert_eq!(reference.context.expect("expected context").id, "NM_004006.2");
}

#[test]
fn parses_nucleotide_substitutions() {
    let genomic = parse_variant("NC_000023.11:g.33038255C>A").into_genomic_edit();
    let cdna = parse_variant("NG_012232.1(NM_004006.2):c.93+1G>T").into_cdna_edit();
    let rna = parse_variant("NM_004006.3:r.76a>u").into_rna_outcome();

    assert_eq!(genomic.location.known_start().coordinate(), Some(33038255));
    assert_nucleotide_substitution(&genomic, "C", "A");

    assert_eq!(cdna.location.known_start().coordinate(), Some(93));
    assert_eq!(cdna.location.known_start().offset(), Some(1));
    assert_nucleotide_substitution(&cdna, "G", "T");

    let (rna_edit, certainty) = rna.produced_edit();
    assert_eq!(*certainty, OutcomeCertainty::Certain);
    assert_eq!(rna_edit.location.known_start().coordinate(), Some(76));
    assert_nucleotide_substitution(rna_edit, "a", "u");
}

#[test]
fn parses_nucleotide_edit_families() {
    let deletion = parse_variant("NM_004006.2:c.5697del").into_cdna_edit();
    let duplication = parse_variant("NC_000001.11:g.1234_2345dup").into_genomic_edit();
    let inversion = parse_variant("NC_000023.10:g.32361330_32361333inv").into_genomic_edit();
    let insertion = parse_variant("LRG_199t1:c.419_420ins[T;450_470;AGGG]").into_cdna_edit();
    let delins = parse_variant("NM_004006.2:c.812_829delinsN[12]").into_cdna_edit();

    assert!(matches!(deletion.kind, NucleotideEditKind::Deletion));
    assert!(matches!(duplication.kind, NucleotideEditKind::Duplication));
    assert!(matches!(inversion.kind, NucleotideEditKind::Inversion));
    assert_eq!(insertion.kind.insertion_items().len(), 3);
    assert!(matches!(
        delins.kind.delins_items()[0],
        NucleotideSequenceItem::Repeat(_)
    ));
}

#[test]
fn parses_cdna_utr_anchored_coordinates() {
    let five_prime_utr = parse_variant("NM_007373.4:c.-81C>T").into_cdna_edit();
    let three_prime_utr = parse_variant("NM_001272071.2:c.*1C>T").into_cdna_edit();
    let five_prime = parse_variant("NM_001385026.1:c.-666+629C>T").into_cdna_edit();
    let three_prime =
        parse_variant("ENSG00000050628.16(ENST00000351052.5):c.*24-12888C>T").into_cdna_edit();

    assert_five_prime_utr_coord(five_prime_utr.location.known_start());
    assert_three_prime_utr_coord(three_prime_utr.location.known_start());

    let five_prime_start = five_prime.location.known_start();
    assert_eq!(
        five_prime_start.anchor(),
        Some(NucleotideAnchor::RelativeCdsStart)
    );
    assert_eq!(five_prime_start.coordinate(), Some(-666));
    assert_eq!(five_prime_start.offset(), Some(629));
    assert_five_prime_intron_coord(five_prime_start);

    let three_prime_start = three_prime.location.known_start();
    assert_eq!(
        three_prime_start.anchor(),
        Some(NucleotideAnchor::RelativeCdsEnd)
    );
    assert_eq!(three_prime_start.coordinate(), Some(24));
    assert_eq!(three_prime_start.offset(), Some(-12888));
    assert_three_prime_intron_coord(three_prime_start);
}

#[test]
fn parses_remote_copied_sequence_in_genomic_insertion() {
    let edit =
        parse_variant("NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]")
            .into_genomic_edit();

    let items = edit.kind.insertion_items();
    assert_eq!(items.len(), 1);
    let NucleotideSequenceItem::Copied(item) = &items[0] else {
        panic!("expected copied sequence item");
    };

    assert_eq!(
        item.source_reference
            .as_ref()
            .expect("expected remote reference")
            .primary
            .id,
        "NC_000022.10"
    );
    assert_eq!(
        item.source_coordinate_system,
        Some(CoordinateSystem::Genomic)
    );
    assert_eq!(item.source_location.start.coordinate(), Some(35788169));
    assert_eq!(item.source_location.interval_end().coordinate(), Some(35788352));
}

#[test]
fn parses_nucleotide_uncertain_locations() {
    let dna = parse_variant("NC_000023.10:g.(33038277_33038278)C>T").into_genomic_edit();
    let rna = parse_variant("NM_004006.2:r.(71_72)_(90_91)del").into_rna_outcome();

    assert!(dna.location.is_uncertain());
    assert_eq!(dna.location.left_bp().start.coordinate(), Some(33038277));
    assert_eq!(
        dna.location.left_bp().interval_end().coordinate(),
        Some(33038278)
    );
    assert_nucleotide_substitution(&dna, "C", "T");

    let (rna_edit, certainty) = rna.produced_edit();
    assert_eq!(*certainty, OutcomeCertainty::Certain);
    assert!(rna_edit.location.is_uncertain());
    assert!(matches!(rna_edit.kind, NucleotideEditKind::Deletion));
    assert_eq!(rna_edit.location.left_bp().start.coordinate(), Some(71));
    assert_eq!(rna_edit.location.right_bp().start.coordinate(), Some(90));
}

#[test]
fn parses_repeat_variants() {
    let genomic = parse_variant("NC_000014.8:g.123_191CAG[19]CAA[4]").into_genomic_edit();
    let rna = parse_variant("NM_004006.3:r.-128_-126[(600_800)]").into_rna_outcome();
    let protein = parse_variant("NP_0123456.1:p.Ala2[10]").into_protein_outcome();

    assert_eq!(
        genomic.kind.repeat_blocks(),
        &[
            make_known_repeat_edit("CAG", 19),
            make_known_repeat_edit("CAA", 4),
        ]
    );

    let (rna_edit, _) = rna.produced_edit();
    assert_eq!(
        rna_edit.kind.repeat_blocks(),
        &[make_shorthand_repeat_edit(make_quantity_range(600, 800))]
    );

    assert_eq!(get_protein_repeat(&protein), &make_shorthand_repeat_edit(10));
}

#[test]
fn parses_repeat_quantity_edges() {
    let unknown = parse_variant("NC_000023.10:g.32717298_32717299insN[?]").into_genomic_edit();
    let lower_unknown =
        parse_variant("NC_000003.12:g.63912687AGC[(?_60)]").into_genomic_edit();
    let upper_unknown =
        parse_variant("NC_000003.12:g.63912687AGC[(60_?)]").into_genomic_edit();

    assert_eq!(
        unknown.kind.insertion_items(),
        &[NucleotideSequenceItem::Repeat(make_unknown_repeat_edit(
            make_unknown_quantity()
        ))]
    );
    assert_eq!(
        lower_unknown.kind.repeat_blocks(),
        &[make_known_repeat_edit_with_quantity(
            "AGC",
            (None, Some(60))
        )]
    );
    assert_eq!(
        upper_unknown.kind.repeat_blocks(),
        &[make_known_repeat_edit_with_quantity(
            "AGC",
            (Some(60), None)
        )]
    );
}

#[test]
fn parses_repeat_sequence_items() {
    let delins = parse_variant("NM_004006.2:c.812_829delinsN[12]").into_cdna_edit();

    assert_eq!(
        delins.kind.delins_items(),
        &[NucleotideSequenceItem::Repeat(make_unknown_repeat_edit(12))]
    );
}

#[test]
fn parses_rna_special_outcomes() {
    let unknown = parse_variant("NM_004006.3:r.?").into_rna_outcome();
    let indeterminate = parse_variant("NM_004006.3:r.(?)").into_rna_outcome();
    let no_change = parse_variant("NM_004006.3:r.=").into_rna_outcome();
    let predicted_no_change = parse_variant("NM_004006.3:r.(=)").into_rna_outcome();
    let none_produced = parse_variant("NM_004006.3:r.0").into_rna_outcome();
    let predicted_none_produced = parse_variant("NM_004006.3:r.0?").into_rna_outcome();
    let splicing = parse_variant("NM_004006.3:r.spl").into_rna_outcome();

    assert!(matches!(unknown, RnaOutcome::Unknown));
    assert!(matches!(indeterminate, RnaOutcome::Indeterminate));
    assert!(matches!(
        no_change,
        RnaOutcome::NoChange(OutcomeCertainty::Certain)
    ));
    assert!(matches!(
        predicted_no_change,
        RnaOutcome::NoChange(OutcomeCertainty::Predicted)
    ));
    assert!(matches!(
        none_produced,
        RnaOutcome::NoneProduced(OutcomeCertainty::Certain)
    ));
    assert!(matches!(
        predicted_none_produced,
        RnaOutcome::NoneProduced(OutcomeCertainty::Predicted)
    ));
    assert!(matches!(splicing, RnaOutcome::UncertainSplicing));
}

#[test]
fn parses_rna_predicted_variant_outcome() {
    let outcome = parse_variant("NM_004006.3:r.(76a>c)").into_rna_outcome();

    let (edit, certainty) = outcome.produced_edit();
    assert_eq!(*certainty, OutcomeCertainty::Predicted);
    assert_eq!(edit.location.known_start().coordinate(), Some(76));
    assert_nucleotide_substitution(edit, "a", "c");
}

#[test]
fn parses_coordinate_specific_alleles() {
    let genomic_form = parse_variant("NC_000001.11:g.[123G>A;345del]").into_genomic_allele_form();
    let cdna_form = parse_variant("NM_004006.2:c.[2376G>C];[2376=]").into_cdna_allele_form();
    let rna_form = parse_variant("NM_004006.3:r.[76a>u];[?]").into_rna_allele_form();

    let AlleleForm::Single(genomic) = genomic_form else {
        panic!("expected a single genomic allele form");
    };
    let AlleleForm::Single(cdna) = cdna_form else {
        panic!("expected a single coding-DNA allele form");
    };
    let AlleleForm::Single(rna) = rna_form else {
        panic!("expected a single RNA allele form");
    };

    assert_eq!(genomic.allele_one.variants.len(), 2);
    assert_eq!(genomic.phase, None);

    assert_eq!(cdna.phase, Some(AllelePhase::Trans));
    assert_eq!(cdna.allele_one.variants.len(), 1);
    assert_eq!(cdna.allele_two.expect("expected second allele").variants.len(), 1);

    assert_eq!(rna.phase, Some(AllelePhase::Trans));
    assert!(matches!(
        rna.allele_two.expect("expected second allele").variants[0],
        RnaOutcome::Unknown
    ));
}

#[test]
fn parses_rna_and_protein_derived_allele_forms() {
    let rna =
        parse_variant("NM_004006.3:r.[897u>g,832_960del,950a>g]").into_rna_allele_form();
    let protein = parse_variant("NP_003997.1:p.[Lys31Asn,Val25_Lys31del,Ser68Arg]")
        .into_protein_allele_form();

    let AlleleForm::Derived(rna) = rna else {
        panic!("expected a derived RNA allele form");
    };
    let AlleleForm::Derived(protein) = protein else {
        panic!("expected a derived protein allele form");
    };

    assert_eq!(rna.outcomes.len(), 3);
    assert!(rna
        .outcomes
        .iter()
        .all(|outcome| matches!(outcome, RnaOutcome::Produced { .. })));

    assert_eq!(protein.outcomes.len(), 3);
    assert!(protein
        .outcomes
        .iter()
        .all(|outcome| matches!(outcome, ProteinOutcome::Produced { .. })));
}

#[test]
fn parses_protein_alternative_allele_form() {
    let protein =
        parse_variant("NP_003997.2:p.[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]")
            .into_protein_allele_form();

    let AlleleForm::Alternative(alternatives) = protein else {
        panic!("expected an alternative protein allele form");
    };

    assert_eq!(alternatives.len(), 2);
    assert_eq!(alternatives[0].phase, Some(AllelePhase::Uncertain));
    assert_eq!(alternatives[1].phase, None);
}

#[test]
fn rejects_malformed_comma_and_alternative_allele_forms() {
    for input in [
        "NP_003997.1:p.[Lys31Asn,]",
        "NP_003997.1:p.[,Lys31Asn]",
        "NP_003997.1:p.[Lys31Asn,,Val25_Lys31del]",
        "NP_003997.1:p.[(Ser68Arg)]^",
        "NP_003997.1:p.^[(Ser68Arg)]",
        "NP_003997.1:p.[(Ser68Arg)]^^[(Asn594del)]",
    ] {
        assert!(parse_hgvs(input).is_err(), "{input} should be rejected");
    }
}

#[test]
fn rejects_unsupported_mixed_protein_allele_forms() {
    for input in [
        "NP_003997.1:p.[Lys31Asn,Val25_Lys31del;Ser68Arg]",
        "NP_003997.1:p.[Lys31Asn;Val25_Lys31del,Ser68Arg]",
        "NP_003997.1:p.[Lys31Asn,Val25_Lys31del]^[(Ser68Arg)]",
        "NP_003997.1:p.[(Ser68Arg)]^[Lys31Asn,Val25_Lys31del]",
    ] {
        assert!(parse_hgvs(input).is_err(), "{input} should be rejected");
    }
}

#[test]
fn parses_uncertain_allele_state() {
    let dna = parse_variant("NC_000001.11:g.123G>A(;)345del").into_genomic_allele();
    let rna = parse_variant("NM_004006.3:r.76a>u(;)(103del)").into_rna_allele();

    assert_eq!(dna.phase, Some(AllelePhase::Uncertain));
    assert_eq!(dna.allele_one.state_certainty, AlleleStateCertainty::Certain);
    assert_eq!(
        dna.allele_two
            .expect("expected second allele")
            .state_certainty,
        AlleleStateCertainty::Certain
    );

    assert_eq!(rna.phase, Some(AllelePhase::Uncertain));
    assert_eq!(
        rna.allele_two
            .expect("expected second allele")
            .state_certainty,
        AlleleStateCertainty::Uncertain
    );
}

#[test]
fn parses_protein_special_outcomes() {
    let unknown = parse_variant("NP_003997.1:p.?").into_protein_outcome();
    let none_produced = parse_variant("NP_003997.1:p.0").into_protein_outcome();
    let predicted_none_produced = parse_variant("NP_003997.1:p.0?").into_protein_outcome();

    assert!(matches!(unknown, ProteinOutcome::Unknown));
    assert!(matches!(
        none_produced,
        ProteinOutcome::NoneProduced(OutcomeCertainty::Certain)
    ));
    assert!(matches!(
        predicted_none_produced,
        ProteinOutcome::NoneProduced(OutcomeCertainty::Predicted)
    ));
}

#[test]
fn parses_protein_substitution_and_prediction() {
    let substitution = parse_variant("NP_003997.1:p.Trp24Ter").into_protein_outcome();
    let predicted = parse_variant("NP_003997.1:p.(Trp24Ter)").into_protein_outcome();

    let (edit, certainty) = substitution.produced_edit();
    assert_eq!(*certainty, OutcomeCertainty::Certain);
    assert_eq!(edit.location.known_start().residue, "Trp");
    assert!(matches!(
        &edit.kind,
        ProteinEditKind::Substitution { to } if to == "Ter"
    ));

    let (_, certainty) = predicted.produced_edit();
    assert_eq!(*certainty, OutcomeCertainty::Predicted);
}

#[test]
fn parses_protein_edit_families() {
    let deletion = parse_variant("NP_003997.2:p.Lys23_Val25del").into_protein_outcome();
    let duplication = parse_variant("NP_003997.1:p.Ser68_Arg70dup").into_protein_outcome();
    let insertion = parse_variant("NP_003997.1:p.Val582_Asn583insAla").into_protein_outcome();
    let delins = parse_variant("NP_003997.1:p.Ser68_Arg70delinsGly").into_protein_outcome();

    assert!(matches!(
        deletion.produced_edit().0.kind,
        ProteinEditKind::Deletion
    ));
    assert!(matches!(
        duplication.produced_edit().0.kind,
        ProteinEditKind::Duplication
    ));
    assert!(matches!(
        insertion.produced_edit().0.kind,
        ProteinEditKind::Insertion { .. }
    ));
    assert!(matches!(
        delins.produced_edit().0.kind,
        ProteinEditKind::DeletionInsertion { .. }
    ));
}

#[test]
fn parses_protein_alleles() {
    let allele = parse_variant("NP_003997.1:p.[Ser68Arg];[Ser68=]").into_protein_allele();

    assert_eq!(allele.phase, Some(AllelePhase::Trans));
    let second = &allele.allele_two.expect("expected second allele").variants[0];
    assert!(matches!(
        second.produced_edit().0.kind,
        ProteinEditKind::NoChange(OutcomeCertainty::Certain)
    ));
}

#[test]
fn parses_protein_frameshift_extension_and_uncertain_location() {
    let frameshift = parse_variant("NP_0123456.1:p.Arg97ProfsTer23").into_protein_outcome();
    let extension = parse_variant("NP_003997.2:p.Ter110GlnextTer17").into_protein_outcome();
    let uncertain_location = parse_variant("NP_003997.1:p.(Ala123_Pro131)Ter").into_protein_outcome();

    let (frameshift_edit, _) = frameshift.produced_edit();
    assert!(matches!(
        &frameshift_edit.kind,
        ProteinEditKind::Frameshift {
            to_residue: Some(residue),
            stop,
        } if residue == "Pro"
            && stop.ordinal == Some(23)
            && stop.kind == ProteinFrameshiftStopKind::Known
    ));

    let (extension_edit, _) = extension.produced_edit();
    assert!(matches!(
        &extension_edit.kind,
        ProteinEditKind::Extension(extension)
            if extension.to_terminal == ProteinExtensionTerminal::C
                && extension.to_residue.as_deref() == Some("Gln")
                && extension.terminal_ordinal == Some(17)
    ));

    let (uncertain_edit, _) = uncertain_location.produced_edit();
    assert!(uncertain_edit.location.is_uncertain());
    assert_eq!(uncertain_edit.location.left_bp().start.residue, "Ala");
    assert_eq!(
        uncertain_edit.location.left_bp().interval_end().residue,
        "Pro"
    );
    assert!(matches!(
        &uncertain_edit.kind,
        ProteinEditKind::Substitution { to } if to == "Ter"
    ));
}

#[test]
fn rejects_malformed_location_shapes() {
    for input in [
        "NC_000023.10:g.(?_?)del",
        "NC_000023.10:g.(?_?)_(?_?)del",
        "NM_004006.2:r.(?_?)del",
        "NM_004006.2:r.(?_?)_(?_?)del",
    ] {
        assert!(parse_hgvs(input).is_err(), "{input} should be rejected");
    }
}
