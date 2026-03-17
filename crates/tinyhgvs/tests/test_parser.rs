use tinyhgvs::{
    parse_hgvs, CoordinateSystem, NucleotideEdit, NucleotidePosition, NucleotidePositionAnchor,
    NucleotideSequenceComponent, NucleotideSequenceSegment, ProteinEdit, ProteinEffect,
    SequenceKind, SequenceSource, VariantDescription,
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
        variant.reference.as_ref().unwrap().primary.kind,
        SequenceKind::RefSeqGeneRegion
    );
    assert_eq!(
        variant
            .reference
            .as_ref()
            .unwrap()
            .context
            .as_ref()
            .unwrap()
            .kind,
        SequenceKind::RefSeqCodingTranscript
    );
    assert_eq!(description.location.start.position, Some(93));
    assert_eq!(description.location.start.offset, 1);
    assert_eq!(
        description.edit,
        NucleotideEdit::Substitution {
            reference: "G".to_string(),
            alternate: "T".to_string(),
        }
    );
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
            assert_eq!(value.location.start.position, Some(5697));
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
            assert_eq!(value.location.start.position, Some(1234));
            assert_eq!(value.location.end.as_ref().unwrap().position, Some(2345));
            assert_eq!(value.edit, NucleotideEdit::Duplication);
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match inversion.description {
        VariantDescription::Nucleotide(value) => {
            assert_eq!(value.location.start.position, Some(32361330));
            assert_eq!(
                value.location.end.as_ref().unwrap().position,
                Some(32361333)
            );
            assert_eq!(value.edit, NucleotideEdit::Inversion);
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }
}

#[test]
fn parses_nucleotide_insertion_sequence_components() {
    let current_reference = parse_variant("LRG_199t1:c.419_420ins[T;450_470;AGGG]");
    let remote_reference =
        parse_variant("NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]");

    match current_reference.description {
        VariantDescription::Nucleotide(value) => {
            let NucleotideEdit::Insertion { sequence } = value.edit else {
                panic!("expected insertion edit");
            };

            assert_eq!(sequence.components.len(), 3);
            assert!(matches!(
                &sequence.components[0],
                NucleotideSequenceComponent::Literal(value) if value == "T"
            ));
            assert!(matches!(
                &sequence.components[1],
                NucleotideSequenceComponent::Segment(NucleotideSequenceSegment {
                    source: SequenceSource::CurrentReference,
                    ..
                })
            ));
            assert!(matches!(
                &sequence.components[2],
                NucleotideSequenceComponent::Literal(value) if value == "AGGG"
            ));
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match remote_reference.description {
        VariantDescription::Nucleotide(value) => {
            let NucleotideEdit::Insertion { sequence } = value.edit else {
                panic!("expected insertion edit");
            };

            assert_eq!(sequence.components.len(), 1);
            match &sequence.components[0] {
                NucleotideSequenceComponent::Segment(NucleotideSequenceSegment {
                    source:
                        SequenceSource::OtherReference {
                            reference,
                            coordinate_system,
                        },
                    location,
                    is_inverted,
                }) => {
                    assert_eq!(reference.primary.kind, SequenceKind::RefSeqChromosome);
                    assert_eq!(*coordinate_system, CoordinateSystem::Genomic);
                    assert_eq!(location.start.position, Some(35788169));
                    assert_eq!(location.end.as_ref().unwrap().position, Some(35788352));
                    assert!(!is_inverted);
                }
                other => panic!("expected remote segment, found {other:?}"),
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
            let NucleotideEdit::DeletionInsertion { sequence } = value.edit else {
                panic!("expected deletion-insertion edit");
            };

            assert!(matches!(
                sequence.components.first().unwrap(),
                NucleotideSequenceComponent::Segment(NucleotideSequenceSegment {
                    source: SequenceSource::CurrentReference,
                    ..
                })
            ));
        }
        other => panic!("expected nucleotide variant, found {other:?}"),
    }

    match repeat.description {
        VariantDescription::Nucleotide(value) => {
            let NucleotideEdit::DeletionInsertion { sequence } = value.edit else {
                panic!("expected deletion-insertion edit");
            };

            assert!(matches!(
                sequence.components.first().unwrap(),
                NucleotideSequenceComponent::Repeat { unit, count }
                    if unit == "N" && *count == 12
            ));
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
                assert_eq!(
                    edit,
                    ProteinEdit::Substitution {
                        to: "Ter".to_string()
                    }
                );
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
}

#[test]
fn normalizes_intronic_and_utr_relative_positions() {
    let intronic = NucleotidePosition {
        anchor: NucleotidePositionAnchor::Coordinate,
        position: Some(93),
        offset: 1,
    };
    let upstream_intronic = NucleotidePosition {
        anchor: NucleotidePositionAnchor::Coordinate,
        position: Some(264),
        offset: -2,
    };
    let five_prime_utr = NucleotidePosition {
        anchor: NucleotidePositionAnchor::CdsStart,
        position: Some(0),
        offset: -81,
    };
    let three_prime_utr = NucleotidePosition {
        anchor: NucleotidePositionAnchor::CdsEnd,
        position: None,
        offset: 1,
    };

    assert!(intronic.is_intronic());
    assert!(!intronic.is_five_prime_utr());
    assert!(!intronic.is_three_prime_utr());

    assert!(upstream_intronic.is_intronic());
    assert!(!five_prime_utr.is_intronic());
    assert!(five_prime_utr.is_five_prime_utr());
    assert!(!five_prime_utr.is_three_prime_utr());

    assert!(!three_prime_utr.is_intronic());
    assert!(!three_prime_utr.is_five_prime_utr());
    assert!(three_prime_utr.is_three_prime_utr());
}

#[test]
fn parses_trimmed_inputs_and_utr_positions() {
    let five_prime = parse_variant("  NM_007373.4:c.-1C>T  ");
    let three_prime = parse_variant("NM_001272071.2:c.*1C>T");
    let upstream_intronic = parse_variant("NG_012232.1(NM_004006.2):c.264-2A>G");

    let VariantDescription::Nucleotide(five_prime) = five_prime.description else {
        panic!("expected nucleotide variant");
    };
    assert_eq!(
        five_prime.location.start.anchor,
        NucleotidePositionAnchor::CdsStart
    );
    assert_eq!(five_prime.location.start.position, Some(0));
    assert_eq!(five_prime.location.start.offset, -1);

    let VariantDescription::Nucleotide(three_prime) = three_prime.description else {
        panic!("expected nucleotide variant");
    };
    assert_eq!(
        three_prime.location.start.anchor,
        NucleotidePositionAnchor::CdsEnd
    );
    assert_eq!(three_prime.location.start.position, None);
    assert_eq!(three_prime.location.start.offset, 1);

    let VariantDescription::Nucleotide(upstream_intronic) = upstream_intronic.description else {
        panic!("expected nucleotide variant");
    };
    assert_eq!(
        upstream_intronic.location.start.anchor,
        NucleotidePositionAnchor::Coordinate
    );
    assert_eq!(upstream_intronic.location.start.position, Some(264));
    assert_eq!(upstream_intronic.location.start.offset, -2);
}

#[test]
fn rejects_examples_deferred_to_future_work() {
    let deferred = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup";
    let error = parse_hgvs(deferred).unwrap_err();

    assert_eq!(error.code(), "unsupported.uncertain_range");
}
