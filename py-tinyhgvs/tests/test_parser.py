import pytest

from tinyhgvs import (
    Accession,
    CoordinateSystem,
    CopiedSequenceItem,
    HgvsVariant,
    Interval,
    LiteralSequenceItem,
    NucleotideAnchor,
    NucleotideCoordinate,
    NucleotideDeletionInsertionEdit,
    NucleotideInsertionEdit,
    NucleotideRepeatBlock,
    NucleotideRepeatEdit,
    NucleotideSequenceOmittedEdit,
    NucleotideSubstitutionEdit,
    NucleotideVariant,
    ParseHgvsErrorKind,
    ProteinCoordinate,
    ProteinDeletionInsertionEdit,
    ProteinEditEffect,
    ProteinInsertionEdit,
    ProteinNoProteinProducedEffect,
    ProteinRepeatEdit,
    ProteinSequence,
    ProteinSequenceOmittedEdit,
    ProteinSubstitutionEdit,
    ProteinUnknownEffect,
    ProteinVariant,
    ReferenceSpec,
    RepeatSequenceItem,
    TinyHGVSError,
    parse_hgvs,
)
from tinyhgvs._tinyhgvs import _roundtrip_variant


def test_parses_nucleotide_substitution_variants():
    variant = parse_hgvs("NG_012232.1(NM_004006.2):c.93+1G>T")

    assert variant.coordinate_system is CoordinateSystem.CODING_DNA
    assert variant.reference is not None
    assert variant.reference.primary.id == "NG_012232.1"
    assert variant.reference.context is not None
    assert variant.reference.context.id == "NM_004006.2"
    assert variant.description.location.start.coordinate == 93
    assert variant.description.location.start.offset == 1
    assert variant.description.edit == NucleotideSubstitutionEdit(
        reference="G",
        alternate="T",
    )

    trimmed = parse_hgvs("  NG_012232.1(NM_004006.2):c.93+1G>T  ")
    assert trimmed.description == variant.description


def test_parses_nucleotide_no_change_and_deletion_variants():
    no_change = parse_hgvs("NM_004006.2:c.123=")
    deletion = parse_hgvs("NM_004006.2:c.5697del")

    assert no_change.description.edit is NucleotideSequenceOmittedEdit.NO_CHANGE
    assert deletion.description.location.start.coordinate == 5697
    assert deletion.description.edit is NucleotideSequenceOmittedEdit.DELETION

    with pytest.raises(TinyHGVSError) as exc_info:
        parse_hgvs("NM_004006.2:c.5697delA")

    assert exc_info.value.code == "invalid.syntax"
    assert exc_info.value.kind is ParseHgvsErrorKind.INVALID_SYNTAX


def test_parses_nucleotide_duplication_and_inversion_variants():
    duplication = parse_hgvs("NC_000001.11:g.1234_2345dup")
    inversion = parse_hgvs("NC_000023.10:g.32361330_32361333inv")

    assert duplication.description.location.start.coordinate == 1234
    assert duplication.description.location.end is not None
    assert duplication.description.location.end.coordinate == 2345
    assert duplication.description.edit is NucleotideSequenceOmittedEdit.DUPLICATION

    assert inversion.description.location.start.coordinate == 32361330
    assert inversion.description.location.end is not None
    assert inversion.description.location.end.coordinate == 32361333
    assert inversion.description.edit is NucleotideSequenceOmittedEdit.INVERSION


def test_parses_nucleotide_insertion_sequence_items():
    current_reference = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")
    remote_reference = parse_hgvs(
        "NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]"
    )

    current_edit = current_reference.description.edit
    assert isinstance(current_edit, NucleotideInsertionEdit)
    assert len(current_edit.items) == 3
    assert current_edit.items[0] == LiteralSequenceItem(value="T")
    assert isinstance(current_edit.items[1], CopiedSequenceItem)
    assert current_edit.items[1].is_from_same_reference is True
    assert current_edit.items[2] == LiteralSequenceItem(value="AGGG")

    remote_edit = remote_reference.description.edit
    assert isinstance(remote_edit, NucleotideInsertionEdit)
    assert len(remote_edit.items) == 1

    remote_item = remote_edit.items[0]
    assert isinstance(remote_item, CopiedSequenceItem)
    assert remote_item.source_reference is not None
    assert remote_item.source_reference.primary.id == "NC_000022.10"
    assert remote_item.source_coordinate_system is CoordinateSystem.GENOMIC
    assert remote_item.source_location.start.coordinate == 35788169
    assert remote_item.source_location.end is not None
    assert remote_item.source_location.end.coordinate == 35788352
    assert remote_item.is_inverted is False
    assert remote_item.is_from_same_reference is False


def test_parses_nucleotide_delins_sequence_forms():
    local_segment = parse_hgvs("NC_000022.10:g.42522624_42522669delins42536337_42536382")
    repeat = parse_hgvs("NM_004006.2:c.812_829delinsN[12]")

    local_edit = local_segment.description.edit
    assert isinstance(local_edit, NucleotideDeletionInsertionEdit)
    assert isinstance(local_edit.items[0], CopiedSequenceItem)
    assert local_edit.items[0].is_from_same_reference is True

    repeat_edit = repeat.description.edit
    assert isinstance(repeat_edit, NucleotideDeletionInsertionEdit)
    assert repeat_edit.items[0] == RepeatSequenceItem(unit="N", count=12)


def test_parses_nucleotide_repeat_variants():
    dna_repeat = parse_hgvs("NC_000014.8:g.123CAG[23]")
    dna_mixed = parse_hgvs("NC_000014.8:g.123_191CAG[19]CAA[4]")
    rna_position_only = parse_hgvs("NM_004006.3:r.-124_-123[14]")
    rna_sequence_given = parse_hgvs("NM_004006.3:r.-110gcu[6]")
    rna_composite = parse_hgvs("NM_004006.3:r.456_465[4]466_489[9]490_499[3]")

    dna_edit = dna_repeat.description.edit
    assert isinstance(dna_edit, NucleotideRepeatEdit)
    assert dna_repeat.description.location.start.coordinate == 123
    assert dna_edit.blocks == (
        NucleotideRepeatBlock(count=23, unit="CAG", location=None),
    )

    mixed_edit = dna_mixed.description.edit
    assert isinstance(mixed_edit, NucleotideRepeatEdit)
    assert dna_mixed.description.location.end is not None
    assert dna_mixed.description.location.end.coordinate == 191
    assert mixed_edit.blocks == (
        NucleotideRepeatBlock(count=19, unit="CAG", location=None),
        NucleotideRepeatBlock(count=4, unit="CAA", location=None),
    )

    position_only_edit = rna_position_only.description.edit
    assert isinstance(position_only_edit, NucleotideRepeatEdit)
    assert rna_position_only.description.location.start.coordinate == -124
    assert rna_position_only.description.location.end is not None
    assert rna_position_only.description.location.end.coordinate == -123
    assert position_only_edit.blocks == (
        NucleotideRepeatBlock(count=14, unit=None, location=None),
    )

    sequence_given_edit = rna_sequence_given.description.edit
    assert isinstance(sequence_given_edit, NucleotideRepeatEdit)
    assert rna_sequence_given.description.location.start.coordinate == -110
    assert rna_sequence_given.description.location.end is None
    assert sequence_given_edit.blocks == (
        NucleotideRepeatBlock(count=6, unit="gcu", location=None),
    )

    composite_edit = rna_composite.description.edit
    assert isinstance(composite_edit, NucleotideRepeatEdit)
    assert rna_composite.description.location.start.coordinate == 456
    assert rna_composite.description.location.end is not None
    assert rna_composite.description.location.end.coordinate == 499
    assert composite_edit.blocks == (
        NucleotideRepeatBlock(count=4, unit=None, location=None),
        NucleotideRepeatBlock(
            count=9,
            unit=None,
            location=Interval(
                start=NucleotideCoordinate(
                    anchor=NucleotideAnchor.ABSOLUTE,
                    coordinate=466,
                ),
                end=NucleotideCoordinate(
                    anchor=NucleotideAnchor.ABSOLUTE,
                    coordinate=489,
                ),
            ),
        ),
        NucleotideRepeatBlock(
            count=3,
            unit=None,
            location=Interval(
                start=NucleotideCoordinate(
                    anchor=NucleotideAnchor.ABSOLUTE,
                    coordinate=490,
                ),
                end=NucleotideCoordinate(
                    anchor=NucleotideAnchor.ABSOLUTE,
                    coordinate=499,
                ),
            ),
        ),
    )


def test_parses_protein_substitution_and_no_change_variants():
    substitution = parse_hgvs("NP_003997.1:p.Trp24Ter")
    no_change = parse_hgvs("NP_003997.1:p.Cys188=")

    assert isinstance(substitution.description, ProteinVariant)
    assert substitution.description.is_predicted is False
    assert isinstance(substitution.description.effect, ProteinEditEffect)
    assert substitution.description.effect.location.start.residue == "Trp"
    assert substitution.description.effect.location.start.ordinal == 24
    assert substitution.description.effect.edit == ProteinSubstitutionEdit(to="Ter")

    assert isinstance(no_change.description.effect, ProteinEditEffect)
    assert no_change.description.effect.edit is ProteinSequenceOmittedEdit.NO_CHANGE


def test_parses_protein_unknown_and_predicted_effects():
    unknown = parse_hgvs("NP_003997.1:p.?")
    predicted = parse_hgvs("LRG_199p1:p.(Met1?)")
    absent = parse_hgvs("LRG_199p1:p.0")

    assert isinstance(unknown.description, ProteinVariant)
    assert unknown.description.effect == ProteinUnknownEffect()
    assert unknown.description.is_predicted is False

    assert isinstance(predicted.description.effect, ProteinEditEffect)
    assert predicted.description.is_predicted is True
    assert predicted.description.effect.location.start.residue == "Met"
    assert predicted.description.effect.location.start.ordinal == 1
    assert predicted.description.effect.edit is ProteinSequenceOmittedEdit.UNKNOWN

    assert absent.description.effect == ProteinNoProteinProducedEffect()


def test_parses_protein_deletion_duplication_insertion_and_delins_variants():
    deletion = parse_hgvs("NP_003997.2:p.Lys23_Val25del")
    duplication = parse_hgvs("NP_003997.2:p.Val7dup")
    insertion = parse_hgvs("p.Lys2_Gly3insGlnSerLys")
    delins = parse_hgvs("p.Cys28delinsTrpVal")

    assert isinstance(deletion.description.effect, ProteinEditEffect)
    assert deletion.description.effect.location.start.residue == "Lys"
    assert deletion.description.effect.location.end is not None
    assert deletion.description.effect.location.end.residue == "Val"
    assert deletion.description.effect.edit is ProteinSequenceOmittedEdit.DELETION

    assert isinstance(duplication.description.effect, ProteinEditEffect)
    assert duplication.description.effect.edit is ProteinSequenceOmittedEdit.DUPLICATION

    assert isinstance(insertion.description.effect, ProteinEditEffect)
    assert insertion.description.effect.edit == ProteinInsertionEdit(
        sequence=ProteinSequence(residues=("Gln", "Ser", "Lys")),
    )

    assert isinstance(delins.description.effect, ProteinEditEffect)
    assert delins.description.effect.edit == ProteinDeletionInsertionEdit(
        sequence=ProteinSequence(residues=("Trp", "Val")),
    )


def test_parses_protein_repeat_variants():
    repeat = parse_hgvs("NP_0123456.1:p.Arg65_Ser67[12]")

    assert isinstance(repeat.description.effect, ProteinEditEffect)
    assert repeat.description.effect.location.start.residue == "Arg"
    assert repeat.description.effect.location.start.ordinal == 65
    assert repeat.description.effect.location.end is not None
    assert repeat.description.effect.location.end.residue == "Ser"
    assert repeat.description.effect.location.end.ordinal == 67
    assert repeat.description.effect.edit == ProteinRepeatEdit(count=12)


def test_distinguishes_intronic_and_utr_relative_coordinates():
    intronic = NucleotideCoordinate(
        anchor=NucleotideAnchor.ABSOLUTE,
        coordinate=93,
        offset=1,
    )
    upstream_intronic = NucleotideCoordinate(
        anchor=NucleotideAnchor.ABSOLUTE,
        coordinate=264,
        offset=-2,
    )
    five_prime_utr = NucleotideCoordinate(
        anchor=NucleotideAnchor.RELATIVE_CDS_START,
        coordinate=-81,
        offset=0,
    )
    three_prime_utr = NucleotideCoordinate(
        anchor=NucleotideAnchor.RELATIVE_CDS_END,
        coordinate=1,
        offset=0,
    )

    assert intronic.is_intronic is True
    assert intronic.is_five_prime_utr is False
    assert intronic.is_three_prime_utr is False

    assert upstream_intronic.is_intronic is True
    assert five_prime_utr.is_intronic is False
    assert five_prime_utr.is_five_prime_utr is True
    assert five_prime_utr.is_three_prime_utr is False

    assert three_prime_utr.is_intronic is False
    assert three_prime_utr.is_five_prime_utr is False
    assert three_prime_utr.is_three_prime_utr is True


def test_parses_utr_and_upstream_intronic_coordinates():
    five_prime = parse_hgvs("NM_007373.4:c.-1C>T")
    three_prime = parse_hgvs("NM_001272071.2:c.*1C>T")
    upstream_intronic = parse_hgvs("NG_012232.1(NM_004006.2):c.264-2A>G")

    assert five_prime.description.location.start.anchor is NucleotideAnchor.RELATIVE_CDS_START
    assert five_prime.description.location.start.coordinate == -1
    assert five_prime.description.location.start.offset == 0

    assert three_prime.description.location.start.anchor is NucleotideAnchor.RELATIVE_CDS_END
    assert three_prime.description.location.start.coordinate == 1
    assert three_prime.description.location.start.offset == 0

    assert upstream_intronic.description.location.start.anchor is NucleotideAnchor.ABSOLUTE
    assert upstream_intronic.description.location.start.coordinate == 264
    assert upstream_intronic.description.location.start.offset == -2


def test_roundtrips_parsed_variant_through_python_to_rust_and_back():
    parsed = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")
    roundtripped = _roundtrip_variant(parsed)

    assert roundtripped == parsed

    utr = parse_hgvs("NM_007373.4:c.-1C>T")
    assert _roundtrip_variant(utr) == utr


def test_roundtrips_manually_constructed_nucleotide_variant():
    variant = HgvsVariant(
        reference=ReferenceSpec(
            primary=Accession(id="LRG_199t1", version=None),
            context=None,
        ),
        coordinate_system=CoordinateSystem.CODING_DNA,
        description=NucleotideVariant(
            location=Interval(
                start=NucleotideCoordinate(
                    anchor=NucleotideAnchor.ABSOLUTE,
                    coordinate=419,
                ),
                end=NucleotideCoordinate(
                    anchor=NucleotideAnchor.ABSOLUTE,
                    coordinate=420,
                ),
            ),
            edit=NucleotideInsertionEdit(
                items=(
                    LiteralSequenceItem(value="T"),
                    CopiedSequenceItem(
                        source_reference=None,
                        source_coordinate_system=None,
                        source_location=Interval(
                            start=NucleotideCoordinate(
                                anchor=NucleotideAnchor.ABSOLUTE,
                                coordinate=450,
                            ),
                            end=NucleotideCoordinate(
                                anchor=NucleotideAnchor.ABSOLUTE,
                                coordinate=470,
                            ),
                        ),
                        is_inverted=False,
                    ),
                    LiteralSequenceItem(value="AGGG"),
                ),
            ),
        ),
    )

    assert _roundtrip_variant(variant) == variant


def test_roundtrips_manually_constructed_nucleotide_repeat_variant():
    variant = HgvsVariant(
        reference=ReferenceSpec(
            primary=Accession(id="NM_004006.3", version=3),
            context=None,
        ),
        coordinate_system=CoordinateSystem.RNA,
        description=NucleotideVariant(
            location=Interval(
                start=NucleotideCoordinate(
                    anchor=NucleotideAnchor.ABSOLUTE,
                    coordinate=456,
                ),
                end=NucleotideCoordinate(
                    anchor=NucleotideAnchor.ABSOLUTE,
                    coordinate=499,
                ),
            ),
            edit=NucleotideRepeatEdit(
                blocks=(
                    NucleotideRepeatBlock(count=4),
                    NucleotideRepeatBlock(
                        count=9,
                        location=Interval(
                            start=NucleotideCoordinate(
                                anchor=NucleotideAnchor.ABSOLUTE,
                                coordinate=466,
                            ),
                            end=NucleotideCoordinate(
                                anchor=NucleotideAnchor.ABSOLUTE,
                                coordinate=489,
                            ),
                        ),
                    ),
                    NucleotideRepeatBlock(
                        count=3,
                        location=Interval(
                            start=NucleotideCoordinate(
                                anchor=NucleotideAnchor.ABSOLUTE,
                                coordinate=490,
                            ),
                            end=NucleotideCoordinate(
                                anchor=NucleotideAnchor.ABSOLUTE,
                                coordinate=499,
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )

    assert _roundtrip_variant(variant) == variant


def test_roundtrips_manually_constructed_protein_variant():
    variant = HgvsVariant(
        reference=ReferenceSpec(
            primary=Accession(id="NP_003997.1", version=1),
            context=None,
        ),
        coordinate_system=CoordinateSystem.PROTEIN,
        description=ProteinVariant(
            is_predicted=False,
            effect=ProteinEditEffect(
                location=Interval(
                    start=ProteinCoordinate(residue="Trp", ordinal=24),
                    end=None,
                ),
                edit=ProteinSubstitutionEdit(to="Ter"),
            ),
        ),
    )

    assert _roundtrip_variant(variant) == variant


def test_roundtrips_manually_constructed_protein_repeat_variant():
    variant = HgvsVariant(
        reference=ReferenceSpec(
            primary=Accession(id="NP_0123456.1", version=1),
            context=None,
        ),
        coordinate_system=CoordinateSystem.PROTEIN,
        description=ProteinVariant(
            is_predicted=False,
            effect=ProteinEditEffect(
                location=Interval(
                    start=ProteinCoordinate(residue="Arg", ordinal=65),
                    end=ProteinCoordinate(residue="Ser", ordinal=67),
                ),
                edit=ProteinRepeatEdit(count=12),
            ),
        ),
    )

    assert _roundtrip_variant(variant) == variant


def test_rejects_examples_deferred_to_future_work():
    deferred = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"

    with pytest.raises(TinyHGVSError) as exc_info:
        parse_hgvs(deferred)

    assert exc_info.value.code == "unsupported.uncertain_range"
    assert exc_info.value.kind is ParseHgvsErrorKind.UNSUPPORTED_SYNTAX


@pytest.mark.parametrize(
    ("example", "code", "kind", "fragment"),
    [
        (
            "NC_000001.11:g.[123G>A;345del]",
            "unsupported.allele",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "[",
        ),
        (
            "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup",
            "unsupported.uncertain_range",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "(",
        ),
        (
            "r.-124_-123[14];[18]",
            "unsupported.allele",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "];[",
        ),
        (
            "NC_000023.11:g.pter_qtersup",
            "unsupported.telomeric_position",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "pter",
        ),
        (
            "NC_000011.10:g.1999904_1999946|gom",
            "unsupported.epigenetic_edit",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "|gom",
        ),
        (
            "NM_001385026.1:c.-666+629C>T",
            "unsupported.cdna_offset_anchor",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "-666+629",
        ),
        (
            "NM_004006.3:r.spl",
            "unsupported.rna_special_state",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "r.spl",
        ),
        (
            "NM_004006.2:r.(222_226)insg",
            "unsupported.rna_uncertain_position",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "r.(...)",
        ),
        (
            "r.-128_-126[(600_800)]",
            "unsupported.uncertain_range",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "[(...)]",
        ),
        (
            "NC_000023.11(NM_004006.2):r.[897u>g,832_960del]",
            "unsupported.rna_splicing_outcome",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "r.[...]",
        ),
        (
            "NM_002354.2:r.-358_555::NM_000251.2:r.212_*279",
            "unsupported.rna_adjoined_transcript",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "::",
        ),
        (
            "p.(Gln576SerfsTer21)",
            "unsupported.protein_frameshift",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "fs",
        ),
        (
            "p.(Ter157Lysext*90)",
            "unsupported.protein_extension",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "ext",
        ),
        (
            "p.(Gln18)[(70_80)]",
            "unsupported.protein_uncertain_consequence",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "[(...)]",
        ),
        (
            "p.Arg78_Gly79insXaa[23]",
            "unsupported.protein_insertion_payload",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "Xaa[...]",
        ),
        (
            "p.(Gly719Ala^Ser)",
            "unsupported.protein_uncertain_consequence",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "^",
        ),
        (
            "not-hgvs",
            "invalid.syntax",
            ParseHgvsErrorKind.INVALID_SYNTAX,
            None,
        ),
    ],
)
def test_raises_tinyhgvs_error_with_structured_diagnostics(
    example: str,
    code: str,
    kind: ParseHgvsErrorKind,
    fragment: str | None,
):
    with pytest.raises(TinyHGVSError) as exc_info:
        parse_hgvs(example)

    error = exc_info.value
    assert isinstance(error, ValueError)
    assert error.code == code
    assert error.kind is kind
    assert error.input == example
    assert error.fragment == fragment
    assert error.message
    assert error.parser_version
    assert f"[{code}]" in str(error)
