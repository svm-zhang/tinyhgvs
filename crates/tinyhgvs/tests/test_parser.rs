use tinyhgvs::{
    parse_hgvs, CoordinateSystem, CopiedSequenceItem, LiteralSequenceItem, NucleotideEdit,
    NucleotideSequenceItem, ProteinEdit, ProteinEffect, ProteinExtensionTerminal,
    ProteinFrameshiftStopKind, RepeatSequenceItem, VariantDescription,
};

fn parse_variant(example: &str) -> tinyhgvs::HgvsVariant {
    parse_hgvs(example).unwrap_or_else(|error| panic!("{example} should parse: {error}"))
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
    assert_eq!(description.location.start.coordinate, 93);
    assert_eq!(description.location.start.offset, 1);
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
        five_prime_intronic.location.start.anchor,
        tinyhgvs::NucleotideAnchor::RelativeCdsStart
    );
    assert_eq!(five_prime_intronic.location.start.coordinate, -666);
    assert_eq!(five_prime_intronic.location.start.offset, 629);

    let VariantDescription::Nucleotide(three_prime_intronic) = three_prime_intronic.description
    else {
        panic!("expected nucleotide variant");
    };
    assert_eq!(
        three_prime_intronic.location.start.anchor,
        tinyhgvs::NucleotideAnchor::RelativeCdsEnd
    );
    assert_eq!(three_prime_intronic.location.start.coordinate, 24);
    assert_eq!(three_prime_intronic.location.start.offset, -12888);

    let VariantDescription::Nucleotide(five_prime_interval) = five_prime_interval.description
    else {
        panic!("expected nucleotide variant");
    };
    let start = five_prime_interval.location.start;
    let end = five_prime_interval
        .location
        .end
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
            assert_eq!(value.location.start.coordinate, 5697);
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
            assert_eq!(value.location.start.coordinate, 1234);
            assert_eq!(value.location.end.as_ref().unwrap().coordinate, 2345);
            assert_eq!(value.edit, NucleotideEdit::Duplication);
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match inversion.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(value.location.start.coordinate, 32361330);
            assert_eq!(value.location.end.as_ref().unwrap().coordinate, 32361333);
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
            assert_eq!(value.location.start.coordinate, 123);
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
            assert_eq!(value.location.start.coordinate, 123);
            assert_eq!(value.location.end.as_ref().unwrap().coordinate, 191);
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
            assert_eq!(value.location.start.coordinate, -124);
            assert_eq!(value.location.end.as_ref().unwrap().coordinate, -123);
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
            assert_eq!(value.location.start.coordinate, -110);
            assert!(value.location.end.is_none());
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
            assert_eq!(value.location.start.coordinate, 456);
            assert_eq!(value.location.end.as_ref().unwrap().coordinate, 499);
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
fn parses_protein_substitution_and_no_change_variants() {
    let substitution = parse_variant("NP_003997.1:p.Trp24Ter");
    let no_change = parse_variant("NP_003997.1:p.Cys188=");

    match substitution.description {
        VariantDescription::Protein(value) => match value.effect {
            ProteinEffect::Edit { location, edit } => {
                assert!(!value.is_predicted);
                assert_eq!(location.start.residue, "Trp");
                assert_eq!(location.start.ordinal, 24);
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
                    assert_eq!(location.start.residue, "Met");
                    assert_eq!(location.start.ordinal, 1);
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
                assert_eq!(location.start.residue, "Lys");
                assert_eq!(location.end.as_ref().unwrap().residue, "Val");
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
                assert_eq!(location.start.residue, "Arg");
                assert_eq!(location.end.as_ref().unwrap().residue, "Ser");
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
                assert_eq!(location.start.residue, "Arg");
                assert_eq!(location.start.ordinal, 97);
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
                assert_eq!(location.start.residue, "Ile");
                assert_eq!(location.start.ordinal, 327);
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
                assert_eq!(location.start.residue, "Met");
                assert_eq!(location.start.ordinal, 1);
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
                assert_eq!(location.start.residue, "Ter");
                assert_eq!(location.start.ordinal, 110);
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
                assert_eq!(location.start.residue, "Ter");
                assert_eq!(location.start.ordinal, 110);
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
                assert_eq!(location.start.residue, "Ter");
                assert_eq!(location.start.ordinal, 327);
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
                assert_eq!(location.start.residue, "Ter");
                assert_eq!(location.start.ordinal, 327);
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
    let intronic = intronic.location.start;
    assert!(intronic.is_intronic());
    assert!(!intronic.is_cds_start_anchored());
    assert!(!intronic.is_cds_end_anchored());
    assert!(!intronic.is_five_prime_utr());
    assert!(!intronic.is_three_prime_utr());

    let VariantDescription::Nucleotide(five_prime_intronic) = five_prime_intronic.description
    else {
        panic!("expected nucleotide variant");
    };
    let five_prime_intronic = five_prime_intronic.location.start;
    assert!(five_prime_intronic.is_intronic());
    assert!(five_prime_intronic.is_cds_start_anchored());
    assert!(!five_prime_intronic.is_cds_end_anchored());
    assert!(!five_prime_intronic.is_five_prime_utr());
    assert!(!five_prime_intronic.is_three_prime_utr());

    let VariantDescription::Nucleotide(five_prime_utr) = five_prime_utr.description else {
        panic!("expected nucleotide variant");
    };
    let five_prime_utr = five_prime_utr.location.start;
    assert!(!five_prime_utr.is_intronic());
    assert!(five_prime_utr.is_cds_start_anchored());
    assert!(!five_prime_utr.is_cds_end_anchored());
    assert!(five_prime_utr.is_five_prime_utr());
    assert!(!five_prime_utr.is_three_prime_utr());

    let VariantDescription::Nucleotide(three_prime_intronic) = three_prime_intronic.description
    else {
        panic!("expected nucleotide variant");
    };
    let three_prime_intronic = three_prime_intronic.location.start;
    assert!(three_prime_intronic.is_intronic());
    assert!(!three_prime_intronic.is_cds_start_anchored());
    assert!(three_prime_intronic.is_cds_end_anchored());
    assert!(!three_prime_intronic.is_five_prime_utr());
    assert!(!three_prime_intronic.is_three_prime_utr());

    let VariantDescription::Nucleotide(three_prime_utr) = three_prime_utr.description else {
        panic!("expected nucleotide variant");
    };
    let three_prime_utr = three_prime_utr.location.start;
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
fn rejects_examples_deferred_to_future_work() {
    let deferred = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup";
    let error = parse_hgvs(deferred).unwrap_err();

    assert_eq!(error.code(), "unsupported.uncertain_range");
}
