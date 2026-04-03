import importlib
import importlib.metadata
import sys

import pytest
import tinyhgvs as tinyhgvs_package
from tinyhgvs import (
    CoordinateSystem,
    CopiedSequenceItem,
    LiteralSequenceItem,
    NucleotideAnchor,
    NucleotideDeletionInsertionEdit,
    NucleotideInsertionEdit,
    NucleotideRepeatEdit,
    NucleotideSequenceOmittedEdit,
    ParseHgvsErrorKind,
    ProteinEditEffect,
    ProteinExtensionTerminal,
    ProteinFrameshiftStopKind,
    ProteinSequenceOmittedEdit,
    TinyHGVSError,
    parse_hgvs,
)


def test_public_package_exports_version_and_core_api():
    assert tinyhgvs_package.__version__
    assert "parse_hgvs" in tinyhgvs_package.__all__
    assert "TinyHGVSError" in tinyhgvs_package.__all__


def test_public_package_falls_back_to_unknown_version_when_metadata_is_missing(
    monkeypatch: pytest.MonkeyPatch,
):
    original_package = sys.modules["tinyhgvs"]

    def raise_not_found(_: str) -> str:
        raise importlib.metadata.PackageNotFoundError

    monkeypatch.setattr(importlib.metadata, "version", raise_not_found)
    sys.modules.pop("tinyhgvs", None)

    try:
        fallback_package = importlib.import_module("tinyhgvs")
        assert fallback_package.__version__ == "0+unknown"
    finally:
        sys.modules.pop("tinyhgvs", None)
        sys.modules["tinyhgvs"] = original_package


def test_parses_nucleotide_substitution_variants():
    variant = parse_hgvs("NG_012232.1(NM_004006.2):c.93+1G>T")

    assert variant.coordinate_system is CoordinateSystem.CODING_DNA
    assert variant.reference is not None
    assert variant.reference.primary.id == "NG_012232.1"
    assert variant.reference.context is not None
    assert variant.reference.context.id == "NM_004006.2"
    assert variant.description.location.start.coordinate == 93
    assert variant.description.location.start.offset == 1
    assert variant.description.edit.reference == "G"
    assert variant.description.edit.alternate == "T"
    assert variant.description.edit.kind == "substitution"

    trimmed = parse_hgvs("  NG_012232.1(NM_004006.2):c.93+1G>T  ")
    assert trimmed.description == variant.description


def test_parses_nucleotide_no_change_and_deletion_variants():
    no_change = parse_hgvs("NM_004006.2:c.123=")
    deletion = parse_hgvs("NM_004006.2:c.5697del")

    assert (
        no_change.description.edit is NucleotideSequenceOmittedEdit.NO_CHANGE
    )
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
    assert (
        duplication.description.edit
        is NucleotideSequenceOmittedEdit.DUPLICATION
    )

    assert inversion.description.location.start.coordinate == 32361330
    assert inversion.description.location.end is not None
    assert inversion.description.location.end.coordinate == 32361333
    assert (
        inversion.description.edit is NucleotideSequenceOmittedEdit.INVERSION
    )


def test_parses_nucleotide_insertion_sequence_items():
    current_reference = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")
    remote_reference = parse_hgvs(
        "NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]"
    )

    current_edit = current_reference.description.edit
    assert isinstance(current_edit, NucleotideInsertionEdit)
    assert len(current_edit.items) == 3
    assert isinstance(current_edit.items[0], LiteralSequenceItem)
    assert getattr(current_edit.items[0], "value", None) == "T"
    assert isinstance(current_edit.items[1], CopiedSequenceItem)
    assert current_edit.items[1].is_from_same_reference is True
    assert getattr(current_edit.items[2], "value", None) == "AGGG"

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
    local_segment = parse_hgvs(
        "NC_000022.10:g.42522624_42522669delins42536337_42536382"
    )
    repeat = parse_hgvs("NM_004006.2:c.812_829delinsN[12]")

    local_edit = local_segment.description.edit
    assert isinstance(local_edit, NucleotideDeletionInsertionEdit)
    assert isinstance(local_edit.items[0], CopiedSequenceItem)
    assert local_edit.items[0].is_from_same_reference is True

    repeat_edit = repeat.description.edit
    assert isinstance(repeat_edit, NucleotideDeletionInsertionEdit)
    assert repeat_edit.items[0].unit == "N"
    assert repeat_edit.items[0].count == 12


def test_parses_nucleotide_repeat_variants():
    dna_repeat = parse_hgvs("NC_000014.8:g.123CAG[23]")
    dna_mixed = parse_hgvs("NC_000014.8:g.123_191CAG[19]CAA[4]")
    rna_position_only = parse_hgvs("NM_004006.3:r.-124_-123[14]")
    rna_sequence_given = parse_hgvs("NM_004006.3:r.-110gcu[6]")
    rna_composite = parse_hgvs("NM_004006.3:r.456_465[4]466_489[9]490_499[3]")

    dna_edit = dna_repeat.description.edit
    assert isinstance(dna_edit, NucleotideRepeatEdit)
    assert dna_repeat.description.location.start.coordinate == 123
    assert len(dna_edit.blocks) == 1
    assert dna_edit.blocks[0].count == 23
    assert dna_edit.blocks[0].unit == "CAG"
    assert dna_edit.blocks[0].location is None

    mixed_edit = dna_mixed.description.edit
    assert isinstance(mixed_edit, NucleotideRepeatEdit)
    assert dna_mixed.description.location.end is not None
    assert dna_mixed.description.location.end.coordinate == 191
    assert len(mixed_edit.blocks) == 2
    assert mixed_edit.blocks[0].count == 19
    assert mixed_edit.blocks[0].unit == "CAG"
    assert mixed_edit.blocks[0].location is None
    assert mixed_edit.blocks[1].count == 4
    assert mixed_edit.blocks[1].unit == "CAA"
    assert mixed_edit.blocks[1].location is None

    position_only_edit = rna_position_only.description.edit
    assert isinstance(position_only_edit, NucleotideRepeatEdit)
    assert rna_position_only.description.location.start.coordinate == -124
    assert rna_position_only.description.location.end is not None
    assert rna_position_only.description.location.end.coordinate == -123
    assert len(position_only_edit.blocks) == 1
    assert position_only_edit.blocks[0].count == 14
    assert position_only_edit.blocks[0].unit is None
    assert position_only_edit.blocks[0].location is None

    sequence_given_edit = rna_sequence_given.description.edit
    assert isinstance(sequence_given_edit, NucleotideRepeatEdit)
    assert rna_sequence_given.description.location.start.coordinate == -110
    assert rna_sequence_given.description.location.end is None
    assert len(sequence_given_edit.blocks) == 1
    assert sequence_given_edit.blocks[0].count == 6
    assert sequence_given_edit.blocks[0].unit == "gcu"
    assert sequence_given_edit.blocks[0].location is None

    composite_edit = rna_composite.description.edit
    assert isinstance(composite_edit, NucleotideRepeatEdit)
    assert rna_composite.description.location.start.coordinate == 456
    assert rna_composite.description.location.end is not None
    assert rna_composite.description.location.end.coordinate == 499
    assert len(composite_edit.blocks) == 3
    assert composite_edit.blocks[0].count == 4
    assert composite_edit.blocks[0].unit is None
    assert composite_edit.blocks[0].location is None
    assert composite_edit.blocks[1].count == 9
    assert composite_edit.blocks[1].unit is None
    assert composite_edit.blocks[1].location is not None
    assert composite_edit.blocks[1].location.start.coordinate == 466
    assert composite_edit.blocks[1].location.end is not None
    assert composite_edit.blocks[1].location.end.coordinate == 489
    assert composite_edit.blocks[2].count == 3
    assert composite_edit.blocks[2].unit is None
    assert composite_edit.blocks[2].location is not None
    assert composite_edit.blocks[2].location.start.coordinate == 490
    assert composite_edit.blocks[2].location.end is not None
    assert composite_edit.blocks[2].location.end.coordinate == 499


def test_parses_protein_substitution_and_no_change_variants():
    substitution = parse_hgvs("NP_003997.1:p.Trp24Ter")
    no_change = parse_hgvs("NP_003997.1:p.Cys188=")

    assert substitution.description.is_predicted is False
    assert isinstance(substitution.description.effect, ProteinEditEffect)
    assert substitution.description.effect.location.start.residue == "Trp"
    assert substitution.description.effect.location.start.ordinal == 24
    assert substitution.description.effect.edit.to == "Ter"
    assert substitution.description.effect.edit.kind == "substitution"

    assert isinstance(no_change.description.effect, ProteinEditEffect)
    assert (
        no_change.description.effect.edit
        is ProteinSequenceOmittedEdit.NO_CHANGE
    )


def test_parses_protein_unknown_and_predicted_effects():
    unknown = parse_hgvs("NP_003997.1:p.?")
    predicted = parse_hgvs("LRG_199p1:p.(Met1?)")
    absent = parse_hgvs("LRG_199p1:p.0")

    assert unknown.description.effect.kind == "unknown"
    assert unknown.description.is_predicted is False

    assert isinstance(predicted.description.effect, ProteinEditEffect)
    assert predicted.description.is_predicted is True
    assert predicted.description.effect.location.start.residue == "Met"
    assert predicted.description.effect.location.start.ordinal == 1
    assert (
        predicted.description.effect.edit is ProteinSequenceOmittedEdit.UNKNOWN
    )

    assert absent.description.effect.kind == "no_protein_produced"


def test_parses_protein_deletion_duplication_insertion_and_delins_variants():
    deletion = parse_hgvs("NP_003997.2:p.Lys23_Val25del")
    duplication = parse_hgvs("NP_003997.2:p.Val7dup")
    insertion = parse_hgvs("p.Lys2_Gly3insGlnSerLys")
    delins = parse_hgvs("p.Cys28delinsTrpVal")

    assert isinstance(deletion.description.effect, ProteinEditEffect)
    assert deletion.description.effect.location.start.residue == "Lys"
    assert deletion.description.effect.location.end is not None
    assert deletion.description.effect.location.end.residue == "Val"
    assert (
        deletion.description.effect.edit is ProteinSequenceOmittedEdit.DELETION
    )

    assert isinstance(duplication.description.effect, ProteinEditEffect)
    assert (
        duplication.description.effect.edit
        is ProteinSequenceOmittedEdit.DUPLICATION
    )

    assert isinstance(insertion.description.effect, ProteinEditEffect)
    assert insertion.description.effect.edit.kind == "insertion"
    assert insertion.description.effect.edit.sequence.residues == (
        "Gln",
        "Ser",
        "Lys",
    )

    assert isinstance(delins.description.effect, ProteinEditEffect)
    assert delins.description.effect.edit.kind == "deletion_insertion"
    assert delins.description.effect.edit.sequence.residues == ("Trp", "Val")


def test_parses_protein_repeat_variants():
    repeat = parse_hgvs("NP_0123456.1:p.Arg65_Ser67[12]")

    assert isinstance(repeat.description.effect, ProteinEditEffect)
    assert repeat.description.effect.location.start.residue == "Arg"
    assert repeat.description.effect.location.start.ordinal == 65
    assert repeat.description.effect.location.end is not None
    assert repeat.description.effect.location.end.residue == "Ser"
    assert repeat.description.effect.location.end.ordinal == 67
    assert repeat.description.effect.edit.kind == "repeat"
    assert repeat.description.effect.edit.count == 12


def test_parses_protein_frameshift_variants():
    short = parse_hgvs("NP_0123456.1:p.Arg97fs")
    long = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
    symbolic_stop = parse_hgvs("NP_0123456.1:p.Arg97Profs*23")
    unknown_stop = parse_hgvs("NP_0123456.1:p.Ile327Argfs*?")
    unknown_stop_ter = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer?")
    predicted = parse_hgvs("p.(Arg97fs)")

    assert isinstance(short.description.effect, ProteinEditEffect)
    assert short.description.effect.location.start.residue == "Arg"
    assert short.description.effect.location.start.ordinal == 97
    assert short.description.effect.edit.kind == "frameshift"
    assert short.description.effect.edit.to_residue is None
    assert short.description.effect.edit.stop.ordinal is None
    assert (
        short.description.effect.edit.stop.kind
        is ProteinFrameshiftStopKind.OMITTED
    )

    assert isinstance(long.description.effect, ProteinEditEffect)
    assert long.description.effect.edit.kind == "frameshift"
    assert long.description.effect.edit.to_residue == "Pro"
    assert long.description.effect.edit.stop.ordinal == 23
    assert (
        long.description.effect.edit.stop.kind
        is ProteinFrameshiftStopKind.KNOWN
    )

    assert isinstance(symbolic_stop.description.effect, ProteinEditEffect)
    assert symbolic_stop.description.effect.edit.kind == "frameshift"
    assert symbolic_stop.description.effect.edit.to_residue == "Pro"
    assert symbolic_stop.description.effect.edit.stop.ordinal == 23
    assert (
        symbolic_stop.description.effect.edit.stop.kind
        is ProteinFrameshiftStopKind.KNOWN
    )

    assert isinstance(unknown_stop.description.effect, ProteinEditEffect)
    assert unknown_stop.description.effect.location.start.residue == "Ile"
    assert unknown_stop.description.effect.location.start.ordinal == 327
    assert unknown_stop.description.effect.edit.kind == "frameshift"
    assert unknown_stop.description.effect.edit.to_residue == "Arg"
    assert unknown_stop.description.effect.edit.stop.ordinal is None
    assert (
        unknown_stop.description.effect.edit.stop.kind
        is ProteinFrameshiftStopKind.UNKNOWN
    )

    assert isinstance(unknown_stop_ter.description.effect, ProteinEditEffect)
    assert unknown_stop_ter.description.effect.edit.kind == "frameshift"
    assert unknown_stop_ter.description.effect.edit.to_residue == "Pro"
    assert unknown_stop_ter.description.effect.edit.stop.ordinal is None
    assert (
        unknown_stop_ter.description.effect.edit.stop.kind
        is ProteinFrameshiftStopKind.UNKNOWN
    )

    assert isinstance(predicted.description.effect, ProteinEditEffect)
    assert predicted.description.is_predicted is True
    assert predicted.description.effect.edit.kind == "frameshift"
    assert predicted.description.effect.edit.to_residue is None
    assert predicted.description.effect.edit.stop.ordinal is None
    assert (
        predicted.description.effect.edit.stop.kind
        is ProteinFrameshiftStopKind.OMITTED
    )


def test_parses_protein_extension_variants():
    n_terminal = parse_hgvs("NP_003997.2:p.Met1ext-5")
    predicted_n_terminal = parse_hgvs("p.(Met1ext-8)")
    c_terminal = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")
    c_terminal_symbolic = parse_hgvs("p.*110Glnext*17")
    unknown_stop = parse_hgvs("p.Ter327ArgextTer?")
    unknown_stop_symbolic = parse_hgvs("p.*327Argext*?")

    assert isinstance(n_terminal.description.effect, ProteinEditEffect)
    assert n_terminal.description.effect.location.start.residue == "Met"
    assert n_terminal.description.effect.location.start.ordinal == 1
    assert n_terminal.description.effect.edit.kind == "extension"
    assert (
        n_terminal.description.effect.edit.to_terminal
        is ProteinExtensionTerminal.N
    )
    assert n_terminal.description.effect.edit.to_residue is None
    assert n_terminal.description.effect.edit.terminal_ordinal == -5

    assert isinstance(predicted_n_terminal.description.effect, ProteinEditEffect)
    assert predicted_n_terminal.description.is_predicted is True
    assert (
        predicted_n_terminal.description.effect.edit.to_terminal
        is ProteinExtensionTerminal.N
    )
    assert predicted_n_terminal.description.effect.edit.terminal_ordinal == -8

    assert isinstance(c_terminal.description.effect, ProteinEditEffect)
    assert c_terminal.description.effect.location.start.residue == "Ter"
    assert c_terminal.description.effect.location.start.ordinal == 110
    assert c_terminal.description.effect.edit.kind == "extension"
    assert (
        c_terminal.description.effect.edit.to_terminal
        is ProteinExtensionTerminal.C
    )
    assert c_terminal.description.effect.edit.to_residue == "Gln"
    assert c_terminal.description.effect.edit.terminal_ordinal == 17

    assert isinstance(c_terminal_symbolic.description.effect, ProteinEditEffect)
    assert c_terminal_symbolic.description.effect.location.start.residue == "Ter"
    assert c_terminal_symbolic.description.effect.location.start.ordinal == 110
    assert (
        c_terminal_symbolic.description.effect.edit.to_terminal
        is ProteinExtensionTerminal.C
    )
    assert c_terminal_symbolic.description.effect.edit.to_residue == "Gln"
    assert c_terminal_symbolic.description.effect.edit.terminal_ordinal == 17

    assert isinstance(unknown_stop.description.effect, ProteinEditEffect)
    assert unknown_stop.description.effect.location.start.residue == "Ter"
    assert unknown_stop.description.effect.location.start.ordinal == 327
    assert (
        unknown_stop.description.effect.edit.to_terminal
        is ProteinExtensionTerminal.C
    )
    assert unknown_stop.description.effect.edit.to_residue == "Arg"
    assert unknown_stop.description.effect.edit.terminal_ordinal is None

    assert isinstance(unknown_stop_symbolic.description.effect, ProteinEditEffect)
    assert (
        unknown_stop_symbolic.description.effect.location.start.residue == "Ter"
    )
    assert unknown_stop_symbolic.description.effect.location.start.ordinal == 327
    assert (
        unknown_stop_symbolic.description.effect.edit.to_terminal
        is ProteinExtensionTerminal.C
    )
    assert unknown_stop_symbolic.description.effect.edit.to_residue == "Arg"
    assert unknown_stop_symbolic.description.effect.edit.terminal_ordinal is None


def test_reports_intronic_and_utr_coordinate_properties_from_parsed_variants():
    intronic = parse_hgvs("NM_004006.2:c.93+1G>T").description.location.start
    upstream_intronic = parse_hgvs(
        "NG_012232.1(NM_004006.2):c.264-2A>G"
    ).description.location.start
    five_prime_utr = parse_hgvs(
        "NM_007373.4:c.-81C>T"
    ).description.location.start
    three_prime_utr = parse_hgvs(
        "NM_001272071.2:c.*1C>T"
    ).description.location.start

    assert intronic.is_intronic is True
    assert intronic.is_five_prime_utr is False
    assert intronic.is_three_prime_utr is False

    assert upstream_intronic.is_intronic is True
    assert upstream_intronic.is_five_prime_utr is False
    assert upstream_intronic.is_three_prime_utr is False

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

    assert (
        five_prime.description.location.start.anchor
        is NucleotideAnchor.RELATIVE_CDS_START
    )
    assert five_prime.description.location.start.coordinate == -1
    assert five_prime.description.location.start.offset == 0

    assert (
        three_prime.description.location.start.anchor
        is NucleotideAnchor.RELATIVE_CDS_END
    )
    assert three_prime.description.location.start.coordinate == 1
    assert three_prime.description.location.start.offset == 0

    assert (
        upstream_intronic.description.location.start.anchor
        is NucleotideAnchor.ABSOLUTE
    )
    assert upstream_intronic.description.location.start.coordinate == 264
    assert upstream_intronic.description.location.start.offset == -2


@pytest.mark.parametrize(
    "example",
    [
        "p.Arg97fsTer23",
        "p.Arg97fs*23",
        "p.Arg97fs*?",
        "p.Arg97Profs",
        "p.Arg97ProfsTer",
        "p.Arg97Profs23",
        "p.Ter97fsTer23",
    ],
)
def test_rejects_malformed_protein_frameshift_variants(example: str):
    with pytest.raises(TinyHGVSError) as exc_info:
        parse_hgvs(example)

    assert exc_info.value.code == "invalid.syntax"
    assert exc_info.value.kind is ParseHgvsErrorKind.INVALID_SYNTAX


@pytest.mark.parametrize(
    "example",
    [
        "p.Met1ext5",
        "p.Met1ext+5",
        "p.Ter110extTer17",
        "p.Ter110Glnext17",
        "p.Ter110GlnextTer",
        "p.Ter110GlnextTer-17",
        "p.Met2ext-5",
    ],
)
def test_rejects_malformed_protein_extension_variants(example: str):
    with pytest.raises(TinyHGVSError) as exc_info:
        parse_hgvs(example)

    assert exc_info.value.code == "invalid.syntax"
    assert exc_info.value.kind is ParseHgvsErrorKind.INVALID_SYNTAX


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
