import pytest

from tinyhgvs import (
    CoordinateSystem,
    CurrentReferenceSource,
    LiteralSequenceComponent,
    NucleotideDeletionEdit,
    NucleotideDeletionInsertionEdit,
    NucleotideDuplicationEdit,
    NucleotideInsertionEdit,
    NucleotideInversionEdit,
    NucleotideNoChangeEdit,
    NucleotidePosition,
    NucleotidePositionAnchor,
    NucleotideSubstitutionEdit,
    OtherReferenceSource,
    ProteinDeletionEdit,
    ProteinDeletionInsertionEdit,
    ProteinDuplicationEdit,
    ProteinEditEffect,
    ProteinInsertionEdit,
    ProteinNoChangeEdit,
    ProteinNoProteinProducedEffect,
    ProteinSequence,
    ProteinSubstitutionEdit,
    ProteinUnknownEdit,
    ProteinUnknownEffect,
    ProteinVariant,
    ParseHgvsErrorKind,
    RepeatSequenceComponent,
    SegmentSequenceComponent,
    SequenceKind,
    TinyHGVSError,
    parse_hgvs,
)


def test_parses_nucleotide_substitution_variants():
    variant = parse_hgvs("NG_012232.1(NM_004006.2):c.93+1G>T")

    assert variant.coordinate_system is CoordinateSystem.CODING_DNA
    assert variant.reference is not None
    assert variant.reference.primary.kind is SequenceKind.REFSEQ_GENE_REGION
    assert variant.reference.context is not None
    assert variant.reference.context.kind is SequenceKind.REFSEQ_CODING_TRANSCRIPT
    assert variant.description.location.start.position == 93
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

    assert no_change.description.edit == NucleotideNoChangeEdit()
    assert deletion.description.location.start.position == 5697
    assert deletion.description.edit == NucleotideDeletionEdit()

    with pytest.raises(TinyHGVSError) as exc_info:
        parse_hgvs("NM_004006.2:c.5697delA")

    assert exc_info.value.code == "invalid.syntax"
    assert exc_info.value.kind is ParseHgvsErrorKind.INVALID_SYNTAX


def test_parses_nucleotide_duplication_and_inversion_variants():
    duplication = parse_hgvs("NC_000001.11:g.1234_2345dup")
    inversion = parse_hgvs("NC_000023.10:g.32361330_32361333inv")

    assert duplication.description.location.start.position == 1234
    assert duplication.description.location.end is not None
    assert duplication.description.location.end.position == 2345
    assert duplication.description.edit == NucleotideDuplicationEdit()

    assert inversion.description.location.start.position == 32361330
    assert inversion.description.location.end is not None
    assert inversion.description.location.end.position == 32361333
    assert inversion.description.edit == NucleotideInversionEdit()


def test_parses_nucleotide_insertion_sequence_components():
    current_reference = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")
    remote_reference = parse_hgvs(
        "NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]"
    )

    current_edit = current_reference.description.edit
    assert isinstance(current_edit, NucleotideInsertionEdit)
    assert len(current_edit.sequence.components) == 3
    assert current_edit.sequence.components[0] == LiteralSequenceComponent(value="T")
    assert isinstance(current_edit.sequence.components[1], SegmentSequenceComponent)
    assert isinstance(
        current_edit.sequence.components[1].segment.source,
        CurrentReferenceSource,
    )
    assert current_edit.sequence.components[2] == LiteralSequenceComponent(value="AGGG")

    remote_edit = remote_reference.description.edit
    assert isinstance(remote_edit, NucleotideInsertionEdit)
    assert len(remote_edit.sequence.components) == 1

    remote_component = remote_edit.sequence.components[0]
    assert isinstance(remote_component, SegmentSequenceComponent)
    assert isinstance(remote_component.segment.source, OtherReferenceSource)
    assert (
        remote_component.segment.source.reference.primary.kind
        is SequenceKind.REFSEQ_CHROMOSOME
    )
    assert remote_component.segment.source.coordinate_system is CoordinateSystem.GENOMIC
    assert remote_component.segment.location.start.position == 35788169
    assert remote_component.segment.location.end is not None
    assert remote_component.segment.location.end.position == 35788352
    assert remote_component.segment.is_inverted is False


def test_parses_nucleotide_delins_sequence_forms():
    local_segment = parse_hgvs("NC_000022.10:g.42522624_42522669delins42536337_42536382")
    repeat = parse_hgvs("NM_004006.2:c.812_829delinsN[12]")

    local_edit = local_segment.description.edit
    assert isinstance(local_edit, NucleotideDeletionInsertionEdit)
    assert isinstance(local_edit.sequence.components[0], SegmentSequenceComponent)
    assert isinstance(
        local_edit.sequence.components[0].segment.source,
        CurrentReferenceSource,
    )

    repeat_edit = repeat.description.edit
    assert isinstance(repeat_edit, NucleotideDeletionInsertionEdit)
    assert repeat_edit.sequence.components[0] == RepeatSequenceComponent(
        unit="N",
        count=12,
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
    assert no_change.description.effect.edit == ProteinNoChangeEdit()


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
    assert predicted.description.effect.edit == ProteinUnknownEdit()

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
    assert deletion.description.effect.edit == ProteinDeletionEdit()

    assert isinstance(duplication.description.effect, ProteinEditEffect)
    assert duplication.description.effect.edit == ProteinDuplicationEdit()

    assert isinstance(insertion.description.effect, ProteinEditEffect)
    assert insertion.description.effect.edit == ProteinInsertionEdit(
        sequence=ProteinSequence(residues=("Gln", "Ser", "Lys")),
    )

    assert isinstance(delins.description.effect, ProteinEditEffect)
    assert delins.description.effect.edit == ProteinDeletionInsertionEdit(
        sequence=ProteinSequence(residues=("Trp", "Val")),
    )


def test_distinguishes_intronic_and_utr_relative_positions():
    intronic = NucleotidePosition(
        anchor=NucleotidePositionAnchor.COORDINATE,
        position=93,
        offset=1,
    )
    upstream_intronic = NucleotidePosition(
        anchor=NucleotidePositionAnchor.COORDINATE,
        position=264,
        offset=-2,
    )
    five_prime_utr = NucleotidePosition(
        anchor=NucleotidePositionAnchor.CDS_START,
        position=0,
        offset=-81,
    )
    three_prime_utr = NucleotidePosition(
        anchor=NucleotidePositionAnchor.CDS_END,
        position=None,
        offset=1,
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


def test_parses_utr_and_upstream_intronic_positions_with_normalized_offsets():
    five_prime = parse_hgvs("NM_007373.4:c.-1C>T")
    three_prime = parse_hgvs("NM_001272071.2:c.*1C>T")
    upstream_intronic = parse_hgvs("NG_012232.1(NM_004006.2):c.264-2A>G")

    assert five_prime.description.location.start.anchor is NucleotidePositionAnchor.CDS_START
    assert five_prime.description.location.start.position == 0
    assert five_prime.description.location.start.offset == -1

    assert three_prime.description.location.start.anchor is NucleotidePositionAnchor.CDS_END
    assert three_prime.description.location.start.position is None
    assert three_prime.description.location.start.offset == 1

    assert upstream_intronic.description.location.start.anchor is (
        NucleotidePositionAnchor.COORDINATE
    )
    assert upstream_intronic.description.location.start.position == 264
    assert upstream_intronic.description.location.start.offset == -2


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
            "NC_000014.8:g.123CAG[23]",
            "unsupported.dna_repeat",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "[23]",
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
            "NM_004006.3:r.9495_9497[4]",
            "unsupported.rna_repeat",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "[4]",
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
            "NP_0123456.1:p.Ala2[10]",
            "unsupported.protein_repeat",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "[10]",
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
