use tinyhgvs::{
    parse_hgvs, AllelePhase, CoordinateSystem, CopiedSequenceItem, LiteralSequenceItem, Location,
    NucleotideEdit, NucleotideSequenceItem, ProteinEdit, ProteinEffect, ProteinExtensionTerminal,
    ProteinFrameshiftStopKind, RepeatSequenceItem, VariantDescription,
};

fn parse_variant(example: &str) -> tinyhgvs::HgvsVariant {
    parse_hgvs(example).unwrap_or_else(|error| panic!("{example} should parse: {error}"))
}

fn known_start<T>(location: &Location<T>) -> &T {
    location.start().expect("expected known location start")
}

fn known_end<T>(location: &Location<T>) -> Option<&T> {
    location.end()
}

#[test]
fn parses_nucleotide_substitution_variants() {
    let variant = parse_variant("NG_012232.1(NM_004006.2):c.93+1G>T");
    let description = match variant.description {
        VariantDescription::Nucleotide(value) => value,
        other => panic!("expected nucleotide variant, found {other:?}"),
    };

    assert_eq!(variant.coordinate_system, CoordinateSystem::CodingDna);
    assert_eq!(
        variant.reference.as_ref().unwrap().primary.id,
        "NG_012232.1"
    );
    assert_eq!(
        variant
            .reference
            .as_ref()
            .unwrap()
            .context
            .as_ref()
            .unwrap()
            .id,
        "NM_004006.2"
    );
    assert_eq!(known_start(&description.location).coordinate, 93);
    assert_eq!(known_start(&description.location).offset, 1);
    let NucleotideEdit::Substitution {
        reference,
        alternate,
    } = description.edit
    else {
        panic!("expected substitution edit");
    };
    assert_eq!(reference, "G");
    assert_eq!(alternate, "T");
}

#[test]
fn parses_cdna_offset_anchor_variants() {
    let five_prime_intronic = parse_variant("NM_001385026.1:c.-666+629C>T");
    let three_prime_intronic =
        parse_variant("ENSG00000050628.16(ENST00000351052.5):c.*24-12888C>T");
    let five_prime_interval = parse_variant("ENST00000440857.1:c.-490-342_-490-341del");

    let VariantDescription::Nucleotide(five_prime_intronic) = five_prime_intronic.description
    else {
        panic!("expected nucleotide variant");
    };
    assert_eq!(
        known_start(&five_prime_intronic.location).anchor,
        tinyhgvs::NucleotideAnchor::RelativeCdsStart
    );
    assert_eq!(known_start(&five_prime_intronic.location).coordinate, -666);
    assert_eq!(known_start(&five_prime_intronic.location).offset, 629);

    let VariantDescription::Nucleotide(three_prime_intronic) = three_prime_intronic.description
    else {
        panic!("expected nucleotide variant");
    };
    assert_eq!(
        known_start(&three_prime_intronic.location).anchor,
        tinyhgvs::NucleotideAnchor::RelativeCdsEnd
    );
    assert_eq!(known_start(&three_prime_intronic.location).coordinate, 24);
    assert_eq!(known_start(&three_prime_intronic.location).offset, -12888);

    let VariantDescription::Nucleotide(five_prime_interval) = five_prime_interval.description
    else {
        panic!("expected nucleotide variant");
    };
    let start = known_start(&five_prime_interval.location);
    let end = five_prime_interval
        .location
        .end()
        .expect("expected interval end");
    assert_eq!(start.anchor, tinyhgvs::NucleotideAnchor::RelativeCdsStart);
    assert_eq!(start.coordinate, -490);
    assert_eq!(start.offset, -342);
    assert_eq!(end.anchor, tinyhgvs::NucleotideAnchor::RelativeCdsStart);
    assert_eq!(end.coordinate, -490);
    assert_eq!(end.offset, -341);
}

#[test]
fn parses_nucleotide_no_change_and_deletion_variants() {
    let no_change = parse_variant("NM_004006.2:c.123=");
    let deletion = parse_variant("NM_004006.2:c.5697del");

    match no_change.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(value.edit, NucleotideEdit::NoChange);
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match deletion.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(known_start(&value.location).coordinate, 5697);
            assert_eq!(value.edit, NucleotideEdit::Deletion);
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    assert!(parse_hgvs("NM_004006.2:c.5697delA").is_err());
}

#[test]
fn parses_nucleotide_duplication_and_inversion_variants() {
    let duplication = parse_variant("NC_000001.11:g.1234_2345dup");
    let inversion = parse_variant("NC_000023.10:g.32361330_32361333inv");

    match duplication.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(known_start(&value.location).coordinate, 1234);
            assert_eq!(known_end(&value.location).unwrap().coordinate, 2345);
            assert_eq!(value.edit, NucleotideEdit::Duplication);
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match inversion.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(known_start(&value.location).coordinate, 32361330);
            assert_eq!(known_end(&value.location).unwrap().coordinate, 32361333);
            assert_eq!(value.edit, NucleotideEdit::Inversion);
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }
}

#[test]
fn parses_nucleotide_insertion_sequence_items() {
    let current_reference = parse_variant("LRG_199t1:c.419_420ins[T;450_470;AGGG]");
    let remote_reference =
        parse_variant("NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]");

    match current_reference.description {
        VariantDescription::Nucleotide(value) => {
            let NucleotideEdit::Insertion { items } = value.edit else {
                panic!("expected insertion edit");
            };

            assert_eq!(items.len(), 3);
            assert!(matches!(
                &items[0],
                NucleotideSequenceItem::Literal(LiteralSequenceItem { value }) if value == "T"
            ));
            assert!(matches!(
                &items[1],
                NucleotideSequenceItem::Copied(CopiedSequenceItem {
                    source_reference: None,
                    source_coordinate_system: None,
                    ..
                })
            ));
            assert!(matches!(
                &items[2],
                NucleotideSequenceItem::Literal(LiteralSequenceItem { value }) if value == "AGGG"
            ));
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match remote_reference.description {
        VariantDescription::Nucleotide(value) => {
            let NucleotideEdit::Insertion { items } = value.edit else {
                panic!("expected insertion edit");
            };

            assert_eq!(items.len(), 1);
            match &items[0] {
                NucleotideSequenceItem::Copied(CopiedSequenceItem {
                    source_reference: Some(reference),
                    source_coordinate_system: Some(coordinate_system),
                    source_location,
                    is_inverted,
                }) => {
                    assert_eq!(reference.primary.id, "NC_000022.10");
                    assert_eq!(*coordinate_system, CoordinateSystem::Genomic);
                    assert_eq!(source_location.start.coordinate, 35788169);
                    assert_eq!(source_location.end.as_ref().unwrap().coordinate, 35788352);
                    assert!(!is_inverted);
                }
                other => panic!("expected remote copied item, found {other:?}"),
            }
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }
}

#[test]
fn parses_nucleotide_delins_sequence_forms() {
    let local_segment = parse_variant("NC_000022.10:g.42522624_42522669delins42536337_42536382");
    let repeat = parse_variant("NM_004006.2:c.812_829delinsN[12]");

    match local_segment.description {
        VariantDescription::Nucleotide(value) => {
            let NucleotideEdit::DeletionInsertion { items } = value.edit else {
                panic!("expected deletion-insertion edit");
            };

            assert!(matches!(
                items.first().unwrap(),
                NucleotideSequenceItem::Copied(CopiedSequenceItem {
                    source_reference: None,
                    source_coordinate_system: None,
                    ..
                })
            ));
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match repeat.description {
        VariantDescription::Nucleotide(value) => {
            let NucleotideEdit::DeletionInsertion { items } = value.edit else {
                panic!("expected deletion-insertion edit");
            };

            assert!(matches!(
                items.first().unwrap(),
                NucleotideSequenceItem::Repeat(RepeatSequenceItem { unit, count })
                    if unit == "N" && *count == 12
            ));
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }
}

#[test]
fn parses_nucleotide_repeat_variants() {
    let dna_repeat = parse_variant("NC_000014.8:g.123CAG[23]");
    let dna_mixed = parse_variant("NC_000014.8:g.123_191CAG[19]CAA[4]");
    let rna_position_only = parse_variant("NM_004006.3:r.-124_-123[14]");
    let rna_sequence_given = parse_variant("NM_004006.3:r.-110gcu[6]");
    let rna_composite = parse_variant("NM_004006.3:r.456_465[4]466_489[9]490_499[3]");

    match dna_repeat.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(known_start(&value.location).coordinate, 123);
            let NucleotideEdit::Repeat { blocks } = value.edit else {
                panic!("expected repeat edit");
            };
            assert_eq!(blocks.len(), 1);
            assert_eq!(blocks[0].count, 23);
            assert_eq!(blocks[0].unit.as_deref(), Some("CAG"));
            assert!(blocks[0].location.is_none());
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match dna_mixed.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(known_start(&value.location).coordinate, 123);
            assert_eq!(known_end(&value.location).unwrap().coordinate, 191);
            let NucleotideEdit::Repeat { blocks } = value.edit else {
                panic!("expected repeat edit");
            };
            assert_eq!(blocks.len(), 2);
            assert_eq!(blocks[0].count, 19);
            assert_eq!(blocks[0].unit.as_deref(), Some("CAG"));
            assert!(blocks[0].location.is_none());
            assert_eq!(blocks[1].count, 4);
            assert_eq!(blocks[1].unit.as_deref(), Some("CAA"));
            assert!(blocks[1].location.is_none());
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match rna_position_only.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(known_start(&value.location).coordinate, -124);
            assert_eq!(known_end(&value.location).unwrap().coordinate, -123);
            let NucleotideEdit::Repeat { blocks } = value.edit else {
                panic!("expected repeat edit");
            };
            assert_eq!(blocks.len(), 1);
            assert_eq!(blocks[0].count, 14);
            assert!(blocks[0].unit.is_none());
            assert!(blocks[0].location.is_none());
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match rna_sequence_given.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(known_start(&value.location).coordinate, -110);
            assert!(known_end(&value.location).is_none());
            let NucleotideEdit::Repeat { blocks } = value.edit else {
                panic!("expected repeat edit");
            };
            assert_eq!(blocks.len(), 1);
            assert_eq!(blocks[0].count, 6);
            assert_eq!(blocks[0].unit.as_deref(), Some("gcu"));
            assert!(blocks[0].location.is_none());
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match rna_composite.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(known_start(&value.location).coordinate, 456);
            assert_eq!(known_end(&value.location).unwrap().coordinate, 499);
            let NucleotideEdit::Repeat { blocks } = value.edit else {
                panic!("expected repeat edit");
            };
            assert_eq!(blocks.len(), 3);
            assert_eq!(blocks[0].count, 4);
            assert_eq!(blocks[0].unit, None);
            assert_eq!(blocks[0].location, None);
            assert_eq!(blocks[1].count, 9);
            assert_eq!(blocks[1].location.as_ref().unwrap().start.coordinate, 466);
            assert_eq!(
                blocks[1]
                    .location
                    .as_ref()
                    .unwrap()
                    .end
                    .as_ref()
                    .unwrap()
                    .coordinate,
                489
            );
            assert_eq!(blocks[2].count, 3);
            assert_eq!(blocks[2].location.as_ref().unwrap().start.coordinate, 490);
            assert_eq!(
                blocks[2]
                    .location
                    .as_ref()
                    .unwrap()
                    .end
                    .as_ref()
                    .unwrap()
                    .coordinate,
                499
            );
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }
}

#[test]
fn parses_nucleotide_allele_variants() {
    let cis = parse_variant("NC_000001.11:g.[123G>A;345del]");
    let trans = parse_variant("NM_004006.3:r.[123c>a];[345del]");
    let uncertain = parse_variant("NC_000001.11:g.123G>A(;)345del");
    let unchanged = parse_variant("NM_004006.2:c.[2376G>C];[2376=]");
    let mixed = parse_variant("NM_004006.2:c.[296T>G;476T>C];[476T>C](;)1083A>C");

    let VariantDescription::NucleotideAllele(cis) = cis.description else {
        panic!("expected nucleotide allele");
    };
    assert_eq!(cis.allele_one.variants.len(), 2);
    assert!(cis.allele_two.is_none());
    assert!(cis.phase.is_none());
    assert!(cis.alleles_unphased.is_empty());
    assert_eq!(cis.iter().count(), 1);
    assert!(matches!(
        cis.allele_one.variants[0].edit,
        NucleotideEdit::Substitution { ref reference, ref alternate }
            if reference == "G" && alternate == "A"
    ));
    assert_eq!(
        known_start(&cis.allele_one.variants[1].location).coordinate,
        345
    );
    assert_eq!(cis.allele_one.variants[1].edit, NucleotideEdit::Deletion);

    let VariantDescription::NucleotideAllele(trans) = trans.description else {
        panic!("expected nucleotide allele");
    };
    assert_eq!(trans.allele_one.variants.len(), 1);
    assert_eq!(trans.phase, Some(AllelePhase::Trans));
    let trans_allele_two = trans.allele_two.as_ref().expect("expected allele two");
    assert_eq!(trans_allele_two.variants.len(), 1);
    assert!(trans.alleles_unphased.is_empty());
    assert_eq!(
        known_start(&trans_allele_two.variants[0].location).coordinate,
        345
    );

    let VariantDescription::NucleotideAllele(uncertain) = uncertain.description else {
        panic!("expected nucleotide allele");
    };
    assert_eq!(uncertain.allele_one.variants.len(), 1);
    assert_eq!(uncertain.phase, Some(AllelePhase::Uncertain));
    let uncertain_allele_two = uncertain.allele_two.as_ref().expect("expected allele two");
    assert!(uncertain.alleles_unphased.is_empty());
    assert_eq!(
        known_start(&uncertain_allele_two.variants[0].location).coordinate,
        345
    );
    assert_eq!(
        uncertain_allele_two.variants[0].edit,
        NucleotideEdit::Deletion
    );

    let VariantDescription::NucleotideAllele(unchanged) = unchanged.description else {
        panic!("expected nucleotide allele");
    };
    assert_eq!(unchanged.allele_one.variants.len(), 1);
    assert_eq!(unchanged.phase, Some(AllelePhase::Trans));
    let unchanged_allele_two = unchanged.allele_two.as_ref().expect("expected allele two");
    assert!(unchanged.alleles_unphased.is_empty());
    assert_eq!(
        known_start(&unchanged_allele_two.variants[0].location).coordinate,
        2376
    );
    assert_eq!(
        unchanged_allele_two.variants[0].edit,
        NucleotideEdit::NoChange
    );

    let VariantDescription::NucleotideAllele(mixed) = mixed.description else {
        panic!("expected nucleotide allele");
    };
    assert_eq!(mixed.allele_one.variants.len(), 2);
    assert_eq!(mixed.phase, Some(AllelePhase::Trans));
    let mixed_allele_two = mixed.allele_two.as_ref().expect("expected allele two");
    assert_eq!(mixed.alleles_unphased.len(), 1);
    assert_eq!(
        known_start(&mixed_allele_two.variants[0].location).coordinate,
        476
    );
    assert_eq!(
        known_start(&mixed.alleles_unphased[0].variants[0].location).coordinate,
        1083
    );
}

#[test]
fn reports_nucleotide_allele_helper_views() {
    let cis = parse_variant("NC_000001.11:g.[123G>A;345del]");
    let trans = parse_variant("NM_004006.3:r.[123c>a];[345del]");
    let uncertain = parse_variant("NC_000001.11:g.123G>A(;)345del");
    let mixed = parse_variant("NM_004006.2:c.[296T>G];[476T>C](;)1083G>C(;)1406del");

    let VariantDescription::NucleotideAllele(cis) = cis.description else {
        panic!("expected nucleotide allele");
    };
    assert!(cis.phased_alleles().is_none());
    assert_eq!(cis.unphased_alleles().len(), 0);
    assert_eq!(cis.allele_one.variants.len(), 2);
    assert!(cis.allele_two.is_none());

    let VariantDescription::NucleotideAllele(trans) = trans.description else {
        panic!("expected nucleotide allele");
    };
    let (allele_one, allele_two) = trans.phased_alleles().expect("expected phased alleles");
    assert_eq!(allele_one.variants.len(), 1);
    assert_eq!(allele_two.variants.len(), 1);
    assert_eq!(trans.unphased_alleles().len(), 0);
    assert_eq!(trans.allele_one.variants.len(), 1);
    assert_eq!(
        trans
            .allele_two
            .as_ref()
            .expect("expected second haplotype")
            .variants
            .len(),
        1
    );

    let VariantDescription::NucleotideAllele(uncertain) = uncertain.description else {
        panic!("expected nucleotide allele");
    };
    assert!(uncertain.phased_alleles().is_none());
    assert_eq!(uncertain.unphased_alleles().len(), 0);
    assert_eq!(uncertain.allele_one.variants.len(), 1);
    assert_eq!(
        known_start(
            &uncertain
                .allele_two
                .as_ref()
                .expect("expected second haplotype")
                .variants[0]
                .location
        )
        .coordinate,
        345
    );

    let VariantDescription::NucleotideAllele(mixed) = mixed.description else {
        panic!("expected nucleotide allele");
    };
    assert!(mixed.phased_alleles().is_some());
    assert_eq!(mixed.unphased_alleles().len(), 2);
    assert_eq!(mixed.allele_one.variants.len(), 1);
    assert_eq!(
        known_start(
            &mixed
                .allele_two
                .as_ref()
                .expect("expected second haplotype")
                .variants[0]
                .location
        )
        .coordinate,
        476
    );
    assert_eq!(
        known_start(&mixed.unphased_alleles()[0].variants[0].location).coordinate,
        1083
    );
    assert_eq!(
        known_start(&mixed.unphased_alleles()[1].variants[0].location).coordinate,
        1406
    );
}

#[test]
fn rejects_malformed_nucleotide_allele_syntax() {
    let cases = [
        "NC_000001.11:g.[123G>A](;)345del",
        "NC_000001.11:g.123G>A(;)[345del]",
        "NC_000001.11:g.[123G>A](;)[345del]",
        "NC_000001.11:g.[123G>A;;345del]",
        "NC_000001.11:g.[123G>A](;)",
        "NC_000001.11:g.[123G>A][345del]",
        "NC_000001.11:g.[123G>A;]",
        "NM_004006.3:r.;[123c>a]",
    ];

    for input in cases {
        let error = parse_hgvs(input).unwrap_err();
        assert_eq!(
            error.code(),
            "invalid.syntax",
            "unexpected code for {input}"
        );
    }
}

#[test]
fn parses_protein_allele_variants() {
    let single = parse_variant("p.[Ser73Arg]");
    let cis = parse_variant("NP_003997.1:p.[Ser68Arg;Asn594del]");
    let trans = parse_variant("NP_003997.1:p.[Ser68Arg];[Ser68=]");
    let uncertain = parse_variant("NP_003997.1:p.(Ser73Arg)(;)(Asn103del)");
    let absent = parse_variant("p.[Ser86Arg];[0]");
    let mixed = parse_variant("p.[Phe233Leu;(Cys690Trp)]");
    let whole_predicted = parse_variant("NP_003997.1:p.[(Ser68Arg;Asn594del)]");
    let range_no_change = parse_variant("p.[Ser68_Arg70dup];[Ser68_Arg70=]");

    let VariantDescription::ProteinAllele(single) = single.description else {
        panic!("expected protein allele");
    };
    assert_eq!(single.allele_one.variants.len(), 1);
    assert!(single.allele_two.is_none());
    assert!(single.phased_alleles().is_none());

    let VariantDescription::ProteinAllele(cis) = cis.description else {
        panic!("expected protein allele");
    };
    assert_eq!(cis.allele_one.variants.len(), 2);
    assert!(cis.allele_two.is_none());
    assert_eq!(cis.iter().count(), 1);
    assert!(cis.phased_alleles().is_none());
    assert_eq!(cis.unphased_alleles().len(), 0);

    let first_cis = &cis.allele_one.variants[0];
    assert!(!first_cis.is_predicted);
    let ProteinEffect::Edit { location, edit } = &first_cis.effect else {
        panic!("expected protein edit");
    };
    assert_eq!(known_start(location).residue, "Ser");
    let ProteinEdit::Substitution { to } = edit else {
        panic!("expected substitution edit");
    };
    assert_eq!(to, "Arg");

    let VariantDescription::ProteinAllele(trans) = trans.description else {
        panic!("expected protein allele");
    };
    assert_eq!(trans.phase, Some(AllelePhase::Trans));
    assert!(trans.phased_alleles().is_some());
    assert!(trans.allele_two.is_some());
    let second_trans = &trans.allele_two.as_ref().unwrap().variants[0];
    let ProteinEffect::Edit { location, edit } = &second_trans.effect else {
        panic!("expected protein edit");
    };
    assert_eq!(known_start(location).residue, "Ser");
    assert!(known_end(location).is_none());
    assert_eq!(*edit, ProteinEdit::NoChange);

    let VariantDescription::ProteinAllele(uncertain) = uncertain.description else {
        panic!("expected protein allele");
    };
    assert_eq!(uncertain.phase, Some(AllelePhase::Uncertain));
    assert!(uncertain.allele_two.is_some());
    assert!(uncertain.phased_alleles().is_none());
    assert_eq!(uncertain.unphased_alleles().len(), 0);
    assert!(uncertain.allele_one.variants[0].is_predicted);
    assert!(uncertain.allele_two.as_ref().unwrap().variants[0].is_predicted);

    let VariantDescription::ProteinAllele(absent) = absent.description else {
        panic!("expected protein allele");
    };
    assert_eq!(absent.phase, Some(AllelePhase::Trans));
    assert_eq!(
        absent.allele_two.as_ref().unwrap().variants[0].effect,
        ProteinEffect::NoProteinProduced
    );

    let VariantDescription::ProteinAllele(mixed) = mixed.description else {
        panic!("expected protein allele");
    };
    assert_eq!(mixed.allele_one.variants.len(), 2);
    assert!(!mixed.allele_one.variants[0].is_predicted);
    assert!(mixed.allele_one.variants[1].is_predicted);

    let VariantDescription::ProteinAllele(whole_predicted) = whole_predicted.description else {
        panic!("expected protein allele");
    };
    assert_eq!(whole_predicted.allele_one.variants.len(), 2);
    assert!(whole_predicted
        .allele_one
        .variants
        .iter()
        .all(|variant| variant.is_predicted));

    let VariantDescription::ProteinAllele(range_no_change) = range_no_change.description else {
        panic!("expected protein allele");
    };
    let second_range = &range_no_change.allele_two.as_ref().unwrap().variants[0];
    let ProteinEffect::Edit { location, edit } = &second_range.effect else {
        panic!("expected protein edit");
    };
    assert_eq!(known_start(location).residue, "Ser");
    assert_eq!(known_end(location).unwrap().residue, "Arg");
    assert_eq!(*edit, ProteinEdit::NoChange);
}

#[test]
fn rejects_malformed_protein_allele_syntax() {
    let cases = [
        "p.([Ser68Arg;Asn594del])",
        "p.([Ser68Arg];[Ser68Arg])",
        "p.[Ser68Arg];[=]",
        "p.[Ser73Arg];[]",
        "p.[Ser68Arg](;)Asn594del",
        "p.[Ser73Arg+p.Asn103del]",
        "p.[Ser73Arg;p.Asn103del]",
    ];

    for input in cases {
        let error = parse_hgvs(input).unwrap_err();
        assert_eq!(
            error.code(),
            "invalid.syntax",
            "unexpected code for {input}"
        );
    }
}

#[test]
fn parses_protein_substitution_and_no_change_variants() {
    let substitution = parse_variant("NP_003997.1:p.Trp24Ter");
    let no_change = parse_variant("NP_003997.1:p.Cys188=");

    match substitution.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert!(!value.is_predicted);
                assert_eq!(known_start(&location).residue, "Trp");
                assert_eq!(known_start(&location).ordinal, 24);
                let ProteinEdit::Substitution { to } = edit else {
                    panic!("expected substitution edit");
                };
                assert_eq!(to, "Ter");
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match no_change.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, .. } => {
                assert_eq!(edit, ProteinEdit::NoChange);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }
}

#[test]
fn parses_protein_unknown_and_predicted_effects() {
    let unknown = parse_variant("NP_003997.1:p.?");
    let predicted = parse_variant("LRG_199p1:p.(Met1?)");
    let absent = parse_variant("LRG_199p1:p.0");

    match unknown.description {
        VariantDescription::Protein(value) => {
            assert_eq!(value.effect, ProteinEffect::Unknown);
            assert!(!value.is_predicted);
        }
        other => panic!("expected protein variant, found {other:?}"),
    }

    match predicted.description {
        VariantDescription::Protein(value) => {
            assert!(value.is_predicted);
            match value.effect {
                ProteinEffect::Edit { location, edit } => {
                    assert_eq!(known_start(&location).residue, "Met");
                    assert_eq!(known_start(&location).ordinal, 1);
                    assert_eq!(edit, ProteinEdit::Unknown);
                }
                other => panic!("expected protein edit, found {other:?}"),
            }
        }
        other => panic!("expected protein variant, found {other:?}"),
    }

    match absent.description {
        VariantDescription::Protein(value) => {
            assert_eq!(value.effect, ProteinEffect::NoProteinProduced);
        }
        other => panic!("expected protein variant, found {other:?}"),
    }
}

#[test]
fn parses_protein_deletion_duplication_insertion_and_delins_variants() {
    let deletion = parse_variant("NP_003997.2:p.Lys23_Val25del");
    let duplication = parse_variant("NP_003997.2:p.Val7dup");
    let insertion = parse_variant("p.Lys2_Gly3insGlnSerLys");
    let delins = parse_variant("p.Cys28delinsTrpVal");
    let repeat = parse_variant("NP_0123456.1:p.Arg65_Ser67[12]");

    match deletion.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, location } => {
                assert_eq!(known_start(&location).residue, "Lys");
                assert_eq!(known_end(&location).unwrap().residue, "Val");
                assert_eq!(edit, ProteinEdit::Deletion);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match duplication.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, .. } => {
                assert_eq!(edit, ProteinEdit::Duplication);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match insertion.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, .. } => {
                let ProteinEdit::Insertion { sequence } = edit else {
                    panic!("expected insertion edit");
                };
                assert_eq!(sequence.residues, vec!["Gln", "Ser", "Lys"]);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match delins.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, .. } => {
                let ProteinEdit::DeletionInsertion { sequence } = edit else {
                    panic!("expected deletion-insertion edit");
                };
                assert_eq!(sequence.residues, vec!["Trp", "Val"]);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match repeat.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, location } => {
                assert_eq!(known_start(&location).residue, "Arg");
                assert_eq!(known_end(&location).unwrap().residue, "Ser");
                let ProteinEdit::Repeat { count } = edit else {
                    panic!("expected repeat edit");
                };
                assert_eq!(count, 12);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }
}

#[test]
fn parses_protein_frameshift_variants() {
    let short = parse_variant("NP_0123456.1:p.Arg97fs");
    let long = parse_variant("NP_0123456.1:p.Arg97ProfsTer23");
    let symbolic_stop = parse_variant("NP_0123456.1:p.Arg97Profs*23");
    let unknown_stop = parse_variant("NP_0123456.1:p.Ile327Argfs*?");
    let unknown_stop_ter = parse_variant("NP_0123456.1:p.Arg97ProfsTer?");
    let predicted = parse_variant("p.(Arg97fs)");

    match short.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert!(!value.is_predicted);
                assert_eq!(known_start(&location).residue, "Arg");
                assert_eq!(known_start(&location).ordinal, 97);
                let ProteinEdit::Frameshift { to_residue, stop } = edit else {
                    panic!("expected frameshift edit");
                };
                assert!(to_residue.is_none());
                assert!(stop.ordinal.is_none());
                assert_eq!(stop.kind, ProteinFrameshiftStopKind::Omitted);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match long.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, .. } => {
                let ProteinEdit::Frameshift { to_residue, stop } = edit else {
                    panic!("expected frameshift edit");
                };
                assert_eq!(to_residue.as_deref(), Some("Pro"));
                assert_eq!(stop.ordinal, Some(23));
                assert_eq!(stop.kind, ProteinFrameshiftStopKind::Known);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match symbolic_stop.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, .. } => {
                let ProteinEdit::Frameshift { to_residue, stop } = edit else {
                    panic!("expected frameshift edit");
                };
                assert_eq!(to_residue.as_deref(), Some("Pro"));
                assert_eq!(stop.ordinal, Some(23));
                assert_eq!(stop.kind, ProteinFrameshiftStopKind::Known);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match unknown_stop.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert_eq!(known_start(&location).residue, "Ile");
                assert_eq!(known_start(&location).ordinal, 327);
                let ProteinEdit::Frameshift { to_residue, stop } = edit else {
                    panic!("expected frameshift edit");
                };
                assert_eq!(to_residue.as_deref(), Some("Arg"));
                assert!(stop.ordinal.is_none());
                assert_eq!(stop.kind, ProteinFrameshiftStopKind::Unknown);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match unknown_stop_ter.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, .. } => {
                let ProteinEdit::Frameshift { to_residue, stop } = edit else {
                    panic!("expected frameshift edit");
                };
                assert_eq!(to_residue.as_deref(), Some("Pro"));
                assert!(stop.ordinal.is_none());
                assert_eq!(stop.kind, ProteinFrameshiftStopKind::Unknown);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match predicted.description {
        VariantDescription::Protein(value) => {
            assert!(value.is_predicted);
            match value.effect {
                ProteinEffect::Edit { edit, .. } => {
                    let ProteinEdit::Frameshift { to_residue, stop } = edit else {
                        panic!("expected frameshift edit");
                    };
                    assert!(to_residue.is_none());
                    assert!(stop.ordinal.is_none());
                    assert_eq!(stop.kind, ProteinFrameshiftStopKind::Omitted);
                }
                other => panic!("expected protein edit, found {other:?}"),
            }
        }
        other => panic!("expected protein variant, found {other:?}"),
    }
}

#[test]
fn parses_protein_extension_variants() {
    let n_terminal = parse_variant("NP_003997.2:p.Met1ext-5");
    let predicted_n_terminal = parse_variant("p.(Met1ext-8)");
    let c_terminal = parse_variant("NP_003997.2:p.Ter110GlnextTer17");
    let c_terminal_symbolic = parse_variant("p.*110Glnext*17");
    let unknown_stop = parse_variant("p.Ter327ArgextTer?");
    let unknown_stop_symbolic = parse_variant("p.*327Argext*?");

    match n_terminal.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert!(!value.is_predicted);
                assert_eq!(known_start(&location).residue, "Met");
                assert_eq!(known_start(&location).ordinal, 1);
                let ProteinEdit::Extension(extension) = edit else {
                    panic!("expected extension edit");
                };
                assert_eq!(extension.to_terminal, ProteinExtensionTerminal::N);
                assert!(extension.to_residue.is_none());
                assert_eq!(extension.terminal_ordinal, Some(-5));
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match predicted_n_terminal.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { edit, .. } => {
                assert!(value.is_predicted);
                let ProteinEdit::Extension(extension) = edit else {
                    panic!("expected extension edit");
                };
                assert_eq!(extension.to_terminal, ProteinExtensionTerminal::N);
                assert_eq!(extension.terminal_ordinal, Some(-8));
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match c_terminal.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert_eq!(known_start(&location).residue, "Ter");
                assert_eq!(known_start(&location).ordinal, 110);
                let ProteinEdit::Extension(extension) = edit else {
                    panic!("expected extension edit");
                };
                assert_eq!(extension.to_terminal, ProteinExtensionTerminal::C);
                assert_eq!(extension.to_residue.as_deref(), Some("Gln"));
                assert_eq!(extension.terminal_ordinal, Some(17));
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match c_terminal_symbolic.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert_eq!(known_start(&location).residue, "Ter");
                assert_eq!(known_start(&location).ordinal, 110);
                let ProteinEdit::Extension(extension) = edit else {
                    panic!("expected extension edit");
                };
                assert_eq!(extension.to_terminal, ProteinExtensionTerminal::C);
                assert_eq!(extension.to_residue.as_deref(), Some("Gln"));
                assert_eq!(extension.terminal_ordinal, Some(17));
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match unknown_stop.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert_eq!(known_start(&location).residue, "Ter");
                assert_eq!(known_start(&location).ordinal, 327);
                let ProteinEdit::Extension(extension) = edit else {
                    panic!("expected extension edit");
                };
                assert_eq!(extension.to_terminal, ProteinExtensionTerminal::C);
                assert_eq!(extension.to_residue.as_deref(), Some("Arg"));
                assert_eq!(extension.terminal_ordinal, None);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }

    match unknown_stop_symbolic.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert_eq!(known_start(&location).residue, "Ter");
                assert_eq!(known_start(&location).ordinal, 327);
                let ProteinEdit::Extension(extension) = edit else {
                    panic!("expected extension edit");
                };
                assert_eq!(extension.to_terminal, ProteinExtensionTerminal::C);
                assert_eq!(extension.to_residue.as_deref(), Some("Arg"));
                assert_eq!(extension.terminal_ordinal, None);
            }
            other => panic!("expected protein edit, found {other:?}"),
        },
        other => panic!("expected protein variant, found {other:?}"),
    }
}

#[test]
fn rejects_malformed_protein_frameshift_variants() {
    let malformed = [
        "p.Arg97fsTer23",
        "p.Arg97fs*23",
        "p.Arg97fs*?",
        "p.Arg97Profs",
        "p.Arg97ProfsTer",
        "p.Arg97Profs23",
        "p.Ter97fsTer23",
    ];

    for example in malformed {
        let error = parse_hgvs(example).unwrap_err();
        assert_eq!(
            error.code(),
            "invalid.syntax",
            "unexpected code for {example}"
        );
    }
}

#[test]
fn rejects_malformed_protein_extension_variants() {
    let malformed = [
        "p.Met1ext5",
        "p.Met1ext+5",
        "p.Ter110extTer17",
        "p.Ter110Glnext17",
        "p.Ter110GlnextTer",
        "p.Ter110GlnextTer-17",
        "p.Met2ext-5",
    ];

    for example in malformed {
        let error = parse_hgvs(example).unwrap_err();
        assert_eq!(
            error.code(),
            "invalid.syntax",
            "unexpected code for {example}"
        );
    }
}

#[test]
fn rejects_malformed_cdna_offset_anchor_variants() {
    let malformed = [
        "NM_001385026.1:c.-106+T>A",
        "NM_001385026.1:c.-106++2T>A",
        "NM_001272071.2:c.*639--1G>A",
        "NM_001272071.2:c.*24-12888_+5del",
        "NM_001385026.1:c.-0+2A>G",
        "NM_001272071.2:c.*0-1G>A",
    ];

    for example in malformed {
        let error = parse_hgvs(example).unwrap_err();
        assert_eq!(
            error.code(),
            "invalid.syntax",
            "unexpected code for {example}"
        );
    }
}

#[test]
fn reports_intronic_and_utr_coordinate_properties_from_parsed_variants() {
    let intronic = parse_variant("NM_004006.2:c.93+1G>T");
    let five_prime_intronic = parse_variant("NM_001385026.1:c.-106+2T>A");
    let five_prime_utr = parse_variant("NM_007373.4:c.-81C>T");
    let three_prime_intronic = parse_variant("NM_001272071.2:c.*639-1G>A");
    let three_prime_utr = parse_variant("NM_001272071.2:c.*1C>T");

    let VariantDescription::Nucleotide(intronic) = intronic.description else {
        panic!("expected nucleotide variant");
    };
    let intronic = known_start(&intronic.location);
    assert!(intronic.is_intronic());
    assert!(!intronic.is_cds_start_anchored());
    assert!(!intronic.is_cds_end_anchored());
    assert!(!intronic.is_five_prime_utr());
    assert!(!intronic.is_three_prime_utr());

    let VariantDescription::Nucleotide(five_prime_intronic) = five_prime_intronic.description
    else {
        panic!("expected nucleotide variant");
    };
    let five_prime_intronic = known_start(&five_prime_intronic.location);
    assert!(five_prime_intronic.is_intronic());
    assert!(five_prime_intronic.is_cds_start_anchored());
    assert!(!five_prime_intronic.is_cds_end_anchored());
    assert!(!five_prime_intronic.is_five_prime_utr());
    assert!(!five_prime_intronic.is_three_prime_utr());

    let VariantDescription::Nucleotide(five_prime_utr) = five_prime_utr.description else {
        panic!("expected nucleotide variant");
    };
    let five_prime_utr = known_start(&five_prime_utr.location);
    assert!(!five_prime_utr.is_intronic());
    assert!(five_prime_utr.is_cds_start_anchored());
    assert!(!five_prime_utr.is_cds_end_anchored());
    assert!(five_prime_utr.is_five_prime_utr());
    assert!(!five_prime_utr.is_three_prime_utr());

    let VariantDescription::Nucleotide(three_prime_intronic) = three_prime_intronic.description
    else {
        panic!("expected nucleotide variant");
    };
    let three_prime_intronic = known_start(&three_prime_intronic.location);
    assert!(three_prime_intronic.is_intronic());
    assert!(!three_prime_intronic.is_cds_start_anchored());
    assert!(three_prime_intronic.is_cds_end_anchored());
    assert!(!three_prime_intronic.is_five_prime_utr());
    assert!(!three_prime_intronic.is_three_prime_utr());

    let VariantDescription::Nucleotide(three_prime_utr) = three_prime_utr.description else {
        panic!("expected nucleotide variant");
    };
    let three_prime_utr = known_start(&three_prime_utr.location);
    assert!(!three_prime_utr.is_intronic());
    assert!(!three_prime_utr.is_cds_start_anchored());
    assert!(three_prime_utr.is_cds_end_anchored());
    assert!(!three_prime_utr.is_five_prime_utr());
    assert!(three_prime_utr.is_three_prime_utr());
}

#[test]
fn trims_input_before_parsing() {
    let trimmed = parse_variant("  NM_007373.4:c.-1C>T  ");
    let plain = parse_variant("NM_007373.4:c.-1C>T");

    assert_eq!(trimmed, plain);
}

#[test]
fn parses_nucleotide_uncertain_locations() {
    let unknown_interval = parse_variant("NC_000023.10:g.?_?del");
    let uncertain_position = parse_variant("NC_000023.10:g.(33038277_33038278)C>T");
    let mixed_uncertain_interval = parse_variant("NC_000023.10:g.(33038277_33038278)_(?_?)del");
    let uncertain_interval = parse_variant("NC_000023.10:g.(?_234567)_(345678_?)del");
    let cdna_uncertain = parse_variant("LRG_199t1:c.(71_72)G>A");
    let rna_unknown_interval = parse_variant("NM_004006.2:r.?_?del");
    let rna_uncertain_interval = parse_variant("NM_004006.2:r.(?_87)del");
    let rna_uncertain_substitution = parse_variant("NM_004006.2:r.(71_72)g>a");

    let VariantDescription::Nucleotide(unknown_interval) = unknown_interval.description else {
        panic!("expected nucleotide variant");
    };
    assert!(!unknown_interval.location.is_uncertain());
    assert!(!unknown_interval.location.is_pos());
    assert!(unknown_interval.location.is_interval());
    assert_eq!(
        unknown_interval.location.start().unwrap().coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );
    assert_eq!(
        unknown_interval.location.end().unwrap().coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );
    assert!(unknown_interval.location.l_interval().is_none());
    assert!(unknown_interval.location.r_interval().is_none());

    let VariantDescription::Nucleotide(uncertain_position) = uncertain_position.description else {
        panic!("expected nucleotide variant");
    };
    assert!(uncertain_position.location.is_uncertain());
    assert!(!uncertain_position.location.is_pos());
    assert!(uncertain_position.location.is_interval());
    assert!(uncertain_position.location.start().is_none());
    assert!(uncertain_position.location.end().is_none());
    let left_region = uncertain_position
        .location
        .l_interval()
        .expect("expected left uncertain region");
    assert_eq!(left_region.start.coordinate, 33038277);
    assert_eq!(left_region.end.as_ref().unwrap().coordinate, 33038278);
    assert!(uncertain_position.location.r_interval().is_none());

    let VariantDescription::Nucleotide(mixed_uncertain_interval) =
        mixed_uncertain_interval.description
    else {
        panic!("expected nucleotide variant");
    };
    assert!(mixed_uncertain_interval.location.is_uncertain());
    let left_region = mixed_uncertain_interval
        .location
        .l_interval()
        .expect("expected left uncertain region");
    assert_eq!(left_region.start.coordinate, 33038277);
    assert_eq!(left_region.end.as_ref().unwrap().coordinate, 33038278);
    let right_region = mixed_uncertain_interval
        .location
        .r_interval()
        .expect("expected right uncertain region");
    assert_eq!(
        right_region.start.coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );
    assert_eq!(
        right_region.end.as_ref().unwrap().coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );

    let VariantDescription::Nucleotide(uncertain_interval) = uncertain_interval.description else {
        panic!("expected nucleotide variant");
    };
    assert!(uncertain_interval.location.is_uncertain());
    assert!(uncertain_interval.location.is_interval());
    let left_region = uncertain_interval
        .location
        .l_interval()
        .expect("expected left uncertain region");
    assert_eq!(
        left_region.start.coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );
    assert_eq!(left_region.end.as_ref().unwrap().coordinate, 234567);
    let right_region = uncertain_interval
        .location
        .r_interval()
        .expect("expected right uncertain region");
    assert_eq!(right_region.start.coordinate, 345678);
    assert_eq!(
        right_region.end.as_ref().unwrap().coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );

    let VariantDescription::Nucleotide(cdna_uncertain) = cdna_uncertain.description else {
        panic!("expected nucleotide variant");
    };
    assert!(cdna_uncertain.location.is_uncertain());
    let left_region = cdna_uncertain
        .location
        .l_interval()
        .expect("expected left uncertain region");
    assert_eq!(left_region.start.coordinate, 71);
    assert_eq!(left_region.end.as_ref().unwrap().coordinate, 72);

    let VariantDescription::Nucleotide(rna_unknown_interval) = rna_unknown_interval.description
    else {
        panic!("expected nucleotide variant");
    };
    assert!(!rna_unknown_interval.location.is_uncertain());
    assert_eq!(
        rna_unknown_interval.location.start().unwrap().coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );
    assert_eq!(
        rna_unknown_interval.location.end().unwrap().coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );

    let VariantDescription::Nucleotide(rna_uncertain_interval) = rna_uncertain_interval.description
    else {
        panic!("expected nucleotide variant");
    };
    assert!(rna_uncertain_interval.location.is_uncertain());
    let left_region = rna_uncertain_interval
        .location
        .l_interval()
        .expect("expected left uncertain region");
    assert_eq!(
        left_region.start.coordinate,
        tinyhgvs::CoordinateKind::Unknown
    );
    assert_eq!(left_region.end.as_ref().unwrap().coordinate, 87);

    let VariantDescription::Nucleotide(rna_uncertain_substitution) =
        rna_uncertain_substitution.description
    else {
        panic!("expected nucleotide variant");
    };
    assert!(rna_uncertain_substitution.location.is_uncertain());
    let left_region = rna_uncertain_substitution
        .location
        .l_interval()
        .expect("expected left uncertain region");
    assert_eq!(left_region.start.coordinate, 71);
    assert_eq!(left_region.end.as_ref().unwrap().coordinate, 72);
}

#[test]
fn parses_protein_uncertain_locations() {
    let uncertain_ter = parse_variant("p.(Ala123_Pro131)Ter");
    let uncertain_frameshift = parse_variant("p.(Ala123_Pro131)fs");

    let VariantDescription::Protein(uncertain_ter) = uncertain_ter.description else {
        panic!("expected protein variant");
    };
    let ProteinEffect::Edit { location, edit } = uncertain_ter.effect else {
        panic!("expected protein edit");
    };
    assert!(location.is_uncertain());
    assert!(!location.is_pos());
    assert!(location.is_interval());
    assert!(matches!(edit, ProteinEdit::Substitution { ref to } if to == "Ter"));
    let left_region = location
        .l_interval()
        .expect("expected uncertain protein region");
    assert_eq!(left_region.start.residue, "Ala");
    assert_eq!(left_region.start.ordinal, 123);
    assert_eq!(left_region.end.as_ref().unwrap().residue, "Pro");
    assert_eq!(left_region.end.as_ref().unwrap().ordinal, 131);
    assert!(location.r_interval().is_none());

    let VariantDescription::Protein(uncertain_frameshift) = uncertain_frameshift.description else {
        panic!("expected protein variant");
    };
    let ProteinEffect::Edit { location, edit } = uncertain_frameshift.effect else {
        panic!("expected protein edit");
    };
    assert!(location.is_uncertain());
    assert!(matches!(edit, ProteinEdit::Frameshift { .. }));
    assert!(location.l_interval().is_some());
}

#[test]
fn rejects_unknown_nucleotide_bounds_with_offsets() {
    for malformed in [
        "NC_000023.10:g.?+1C>T",
        "NC_000023.10:g.(?+1_5)del",
        "NC_000023.10:g.(?-2_5)del",
    ] {
        let error = parse_hgvs(malformed).unwrap_err();
        assert_eq!(error.code(), "invalid.syntax");
    }
}

#[test]
fn rejects_malformed_uncertain_range_variants() {
    for malformed in [
        "NC_000023.10:g.(?_?)del",
        "NC_000023.10:g.(?_?)_(?_?)del",
        "NM_004006.2:r.(?_?)del",
        "NM_004006.2:r.(?_?)_(?_?)del",
    ] {
        let error = parse_hgvs(malformed).unwrap_err();
        assert_eq!(error.code(), "invalid.syntax");
    }
}

#[test]
fn still_rejects_examples_deferred_to_future_work() {
    let deferred = "r.-128_-126[(600_800)]";
    let error = parse_hgvs(deferred).unwrap_err();

    assert_eq!(error.code(), "unsupported.uncertain_size");
}
