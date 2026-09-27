import importlib
import importlib.metadata
import sys

import pytest
import tinyhgvs as tinyhgvs_package
from tinyhgvs import (
    AllelePhase,
    AlleleVariant,
    AlternativeAlleleForm,
    CoordinateSystem,
    CopiedSequence,
    DerivedAlleleForm,
    KnownLocation,
    KnownProteinFrameshiftStop,
    KnownProteinInsertion,
    KnownQuantity,
    KnownRepeatUnit,
    LiteralSequence,
    NucleotideDeletion,
    NucleotideDeletionInsertion,
    NucleotideAnchor,
    NucleotideCoordinate,
    NucleotideDuplication,
    NucleotideInsertion,
    NucleotideInversion,
    NucleotideNoChange,
    NucleotideRepeat,
    OmittedProteinFrameshiftStop,
    OutcomeCertainty,
    ParseHgvsErrorKind,
    ProteinExtensionTerminal,
    ProteinDeletion,
    ProteinDeletionInsertion,
    ProteinDuplication,
    ProteinExtension,
    ProteinFrameshift,
    ProteinNoChange,
    ProteinNotProduced,
    ProteinProduced,
    ProteinProducedAlternatives,
    ProteinRepeat,
    ProteinSubstitution,
    ProteinUnknown,
    Repeat,
    RnaIndeterminate,
    RnaNoChange,
    RnaNotProduced,
    RnaProduced,
    RnaUncertainSplicing,
    RnaUnknown,
    TerminatingProteinInsertion,
    UncertainLocation,
    UncertainQuantity,
    UnknownProteinFrameshiftStop,
    UnknownProteinInsertion,
    UnknownQuantity,
    UnknownRepeatUnit,
    TinyHGVSError,
    parse_hgvs,
)


def test_public_package_exports_version_and_core_api():
    assert tinyhgvs_package.__version__
    assert "parse_hgvs" in tinyhgvs_package.__all__
    assert "TinyHGVSError" in tinyhgvs_package.__all__
    assert "Location" in tinyhgvs_package.__all__
    assert "NucleotideCoordinate" in tinyhgvs_package.__all__


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
    assert variant.description.location.start.is_known is True
    assert variant.description.location.start.coordinate == 93
    assert variant.description.location.start.offset == 1
    assert variant.description.reference == "G"
    assert variant.description.alternate == "T"

    trimmed = parse_hgvs("  NG_012232.1(NM_004006.2):c.93+1G>T  ")
    assert trimmed.description == variant.description


def test_parses_cdna_offset_anchor_variants():
    five_prime_intronic = parse_hgvs("NM_001385026.1:c.-666+629C>T")
    three_prime_intronic = parse_hgvs(
        "ENSG00000050628.16(ENST00000351052.5):c.*24-12888C>T"
    )
    five_prime_interval = parse_hgvs(
        "ENST00000440857.1:c.-490-342_-490-341del"
    )

    assert (
        five_prime_intronic.description.location.start.anchor
        is NucleotideAnchor.RELATIVE_CDS_START
    )
    assert five_prime_intronic.description.location.start.coordinate == -666
    assert five_prime_intronic.description.location.start.offset == 629

    assert (
        three_prime_intronic.description.location.start.anchor
        is NucleotideAnchor.RELATIVE_CDS_END
    )
    assert three_prime_intronic.description.location.start.coordinate == 24
    assert three_prime_intronic.description.location.start.offset == -12888

    assert (
        five_prime_interval.description.location.start.anchor
        is NucleotideAnchor.RELATIVE_CDS_START
    )
    assert five_prime_interval.description.location.start.coordinate == -490
    assert five_prime_interval.description.location.start.offset == -342
    assert five_prime_interval.description.location.end is not None
    assert (
        five_prime_interval.description.location.end.anchor
        is NucleotideAnchor.RELATIVE_CDS_START
    )
    assert five_prime_interval.description.location.end.coordinate == -490
    assert five_prime_interval.description.location.end.offset == -341


def test_parses_nucleotide_no_change_and_deletion_variants():
    no_change = parse_hgvs("NM_004006.2:c.123=")
    deletion = parse_hgvs("NM_004006.2:c.5697del")

    assert isinstance(no_change.description, NucleotideNoChange)
    assert deletion.description.location.start.coordinate == 5697
    assert isinstance(deletion.description, NucleotideDeletion)

    with pytest.raises(TinyHGVSError) as exc_info:
        parse_hgvs("NM_004006.2:c.5697delA")

    assert exc_info.value.code == "invalid.syntax"
    assert exc_info.value.kind is ParseHgvsErrorKind.INVALID_SYNTAX


def test_parses_nucleotide_duplication_and_inversion_variants():
    duplication = parse_hgvs("NC_000001.11:g.1234_2345dup")
    inversion = parse_hgvs("NC_000023.10:g.32361330_32361333inv")

    assert isinstance(duplication.description, NucleotideDuplication)
    assert duplication.description.location.start.coordinate == 1234
    assert duplication.description.location.end is not None
    assert duplication.description.location.end.coordinate == 2345

    assert isinstance(inversion.description, NucleotideInversion)
    assert inversion.description.location.start.coordinate == 32361330
    assert inversion.description.location.end is not None
    assert inversion.description.location.end.coordinate == 32361333


def test_reports_known_nucleotide_location_helper_views():
    single = parse_hgvs("NM_004006.2:c.5697del")
    location = single.description.location

    assert isinstance(location, KnownLocation)
    assert location.is_position is True
    assert location.is_interval is False
    assert location.start.is_known is True
    assert location.start.coordinate == 5697
    assert location.end is None

    interval = parse_hgvs("NM_004006.2:c.93_94del")
    location = interval.description.location

    assert isinstance(location, KnownLocation)
    assert location.is_position is False
    assert location.is_interval is True
    assert location.start.is_known is True
    assert location.start.coordinate == 93
    assert location.end is not None
    assert location.end.is_known is True
    assert location.end.coordinate == 94


def test_parses_uncertain_nucleotide_locations():
    unknown_range = parse_hgvs("NC_000023.10:g.?_?del")
    assert isinstance(unknown_range.description.location, KnownLocation)
    assert unknown_range.description.location.is_interval is True
    assert unknown_range.description.location.is_position is False
    assert unknown_range.description.location.start.is_unknown is True
    assert unknown_range.description.location.start.is_known is False
    assert unknown_range.description.location.start.anchor is None
    assert unknown_range.description.location.start.coordinate is None
    assert unknown_range.description.location.start.offset is None
    assert unknown_range.description.location.end is not None
    assert unknown_range.description.location.end.is_unknown is True
    assert unknown_range.description.location.end.anchor is None
    assert unknown_range.description.location.end.coordinate is None
    assert unknown_range.description.location.end.offset is None

    single_region = parse_hgvs("NC_000023.10:g.(33038277_33038278)C>T")
    location = single_region.description.location
    assert isinstance(location, UncertainLocation)
    assert location.is_interval is True
    assert location.is_position is False
    assert location.start.start.is_known is True
    assert location.start.start.coordinate == 33038277
    assert location.start.end is not None
    assert location.start.end.coordinate == 33038278
    assert location.end is None

    mixed_unknown = parse_hgvs("NC_000023.10:g.(?_32238146)_(32984039_?)del")
    location = mixed_unknown.description.location
    assert isinstance(location, UncertainLocation)
    assert location.start.start.is_unknown is True
    assert location.start.start.anchor is None
    assert location.start.start.coordinate is None
    assert location.start.start.offset is None
    assert location.start.end is not None
    assert location.start.end.is_known is True
    assert location.start.end.coordinate == 32238146
    assert location.end is not None
    assert location.end.start.is_known is True
    assert location.end.start.coordinate == 32984039
    assert location.end.end is not None
    assert location.end.end.is_unknown is True
    assert location.end.end.anchor is None
    assert location.end.end.coordinate is None
    assert location.end.end.offset is None

    uncertain_range = parse_hgvs("NM_004006.2:r.(71_72)_(90_91)del")
    location = uncertain_range.description.edit.location
    assert isinstance(location, UncertainLocation)
    assert location.start.start.coordinate == 71
    assert location.start.end is not None
    assert location.start.end.coordinate == 72
    assert location.end is not None
    assert location.end.start.coordinate == 90
    assert location.end.end is not None
    assert location.end.end.coordinate == 91

    rna_insertion = parse_hgvs("NM_004006.2:r.(222_226)insg")
    edit = rna_insertion.description.edit
    location = edit.location
    assert isinstance(location, UncertainLocation)
    assert location.start.start.coordinate == 222
    assert location.start.end is not None
    assert location.start.end.coordinate == 226
    assert isinstance(edit, NucleotideInsertion)
    assert isinstance(edit.sequence[0], LiteralSequence)
    assert edit.sequence[0].value == "g"


def test_rejects_malformed_uncertain_nucleotide_locations():
    cases = [
        "NC_000023.10:g.(?_?)del",
        "NC_000023.10:g.(?_?)_(?_?)del",
        "NM_004006.2:r.(?_?)del",
    ]

    for input_value in cases:
        with pytest.raises(TinyHGVSError) as exc_info:
            parse_hgvs(input_value)

        assert exc_info.value.code == "invalid.syntax"


def test_parses_nucleotide_insertion_sequence_items():
    current_reference = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")
    remote_reference = parse_hgvs(
        "NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]"
    )

    current_edit = current_reference.description
    assert isinstance(current_edit, NucleotideInsertion)
    assert len(current_edit.sequence) == 3
    assert isinstance(current_edit.sequence[0], LiteralSequence)
    assert getattr(current_edit.sequence[0], "value", None) == "T"
    assert isinstance(current_edit.sequence[1], CopiedSequence)
    assert current_edit.sequence[1].is_from_same_reference is True
    assert getattr(current_edit.sequence[2], "value", None) == "AGGG"

    remote_edit = remote_reference.description
    assert isinstance(remote_edit, NucleotideInsertion)
    assert len(remote_edit.sequence) == 1

    remote_item = remote_edit.sequence[0]
    assert isinstance(remote_item, CopiedSequence)
    assert remote_item.reference is not None
    assert remote_item.reference.primary.id == "NC_000022.10"
    assert remote_item.coordinate_system is CoordinateSystem.GENOMIC
    assert remote_item.location.start.coordinate == 35788169
    assert remote_item.location.end is not None
    assert remote_item.location.end.coordinate == 35788352
    assert remote_item.is_inverted is False
    assert remote_item.is_from_same_reference is False


def test_parses_nucleotide_delins_sequence_forms():
    local_segment = parse_hgvs(
        "NC_000022.10:g.42522624_42522669delins42536337_42536382"
    )
    repeat = parse_hgvs("NM_004006.2:c.812_829delinsN[12]")

    local_edit = local_segment.description
    assert isinstance(local_edit, NucleotideDeletionInsertion)
    assert isinstance(local_edit.sequence[0], CopiedSequence)
    assert local_edit.sequence[0].is_from_same_reference is True

    repeat_edit = repeat.description
    assert isinstance(repeat_edit, NucleotideDeletionInsertion)
    assert isinstance(repeat_edit.sequence[0], Repeat)
    assert isinstance(repeat_edit.sequence[0].unit, UnknownRepeatUnit)
    assert repeat_edit.sequence[0].quantity.count == 12


def test_parses_nucleotide_repeat_variants():
    dna_repeat = parse_hgvs("NC_000014.8:g.123CAG[23]")
    dna_mixed = parse_hgvs("NC_000014.8:g.123_191CAG[19]CAA[4]")
    rna_position_only = parse_hgvs("NM_004006.3:r.-124_-123[14]")
    rna_sequence_given = parse_hgvs("NM_004006.3:r.-110gcu[6]")
    rna_composite = parse_hgvs("NM_004006.3:r.456_499us[4]cag[9]gccag[3]")

    dna_edit = dna_repeat.description
    assert isinstance(dna_edit, NucleotideRepeat)
    assert dna_repeat.description.location.start.coordinate == 123
    assert len(dna_edit.sequence) == 1
    assert dna_edit.sequence[0].quantity.count == 23
    assert isinstance(dna_edit.sequence[0].unit, KnownRepeatUnit)
    assert dna_edit.sequence[0].unit.value == "CAG"

    mixed_edit = dna_mixed.description
    assert isinstance(mixed_edit, NucleotideRepeat)
    assert dna_mixed.description.location.end is not None
    assert dna_mixed.description.location.end.coordinate == 191
    assert len(mixed_edit.sequence) == 2
    assert mixed_edit.sequence[0].quantity.count == 19
    assert isinstance(mixed_edit.sequence[0].unit, KnownRepeatUnit)
    assert mixed_edit.sequence[0].unit.value == "CAG"
    assert mixed_edit.sequence[1].quantity.count == 4
    assert isinstance(mixed_edit.sequence[1].unit, KnownRepeatUnit)
    assert mixed_edit.sequence[1].unit.value == "CAA"

    position_only_edit = rna_position_only.description.edit
    assert isinstance(position_only_edit, NucleotideRepeat)
    assert position_only_edit.location.start.coordinate == -124
    assert position_only_edit.location.end is not None
    assert position_only_edit.location.end.coordinate == -123
    assert len(position_only_edit.sequence) == 1
    assert position_only_edit.sequence[0].quantity.count == 14
    assert position_only_edit.sequence[0].unit is None

    sequence_given_edit = rna_sequence_given.description.edit
    assert isinstance(sequence_given_edit, NucleotideRepeat)
    assert sequence_given_edit.location.start.coordinate == -110
    assert sequence_given_edit.location.end is None
    assert len(sequence_given_edit.sequence) == 1
    assert sequence_given_edit.sequence[0].quantity.count == 6
    assert isinstance(sequence_given_edit.sequence[0].unit, KnownRepeatUnit)
    assert sequence_given_edit.sequence[0].unit.value == "gcu"

    composite_edit = rna_composite.description.edit
    assert isinstance(composite_edit, NucleotideRepeat)
    assert composite_edit.location.start.coordinate == 456
    assert composite_edit.location.end is not None
    assert composite_edit.location.end.coordinate == 499
    assert len(composite_edit.sequence) == 3
    assert composite_edit.sequence[0].quantity.count == 4
    assert isinstance(composite_edit.sequence[0].unit, KnownRepeatUnit)
    assert composite_edit.sequence[0].unit.value == "us"
    assert composite_edit.sequence[1].quantity.count == 9
    assert isinstance(composite_edit.sequence[1].unit, KnownRepeatUnit)
    assert composite_edit.sequence[1].unit.value == "cag"
    assert composite_edit.sequence[2].quantity.count == 3
    assert isinstance(composite_edit.sequence[2].unit, KnownRepeatUnit)
    assert composite_edit.sequence[2].unit.value == "gccag"


def test_reports_repeat_quantity_helper_views():
    known = parse_hgvs("NC_000014.8:g.123CAG[23]")
    unknown = parse_hgvs("NC_000003.12:g.63912687AGC[?]")
    uncertain = parse_hgvs("NC_000003.12:g.63912687AGC[(19_23)]")
    omitted_unit = parse_hgvs("NM_004006.3:r.-124_-123[14]")

    repeat = known.description.sequence[0]
    assert isinstance(repeat, Repeat)
    assert repeat.is_unit_known is True
    assert repeat.is_quantity_known is True
    assert repeat.is_quantity_unknown is False
    assert repeat.is_quantity_uncertain is False
    assert isinstance(repeat.quantity, KnownQuantity)
    assert repeat.quantity.count == 23

    repeat = unknown.description.sequence[0]
    assert repeat.is_quantity_unknown is True
    assert repeat.is_quantity_known is False
    assert repeat.is_quantity_uncertain is False
    assert isinstance(repeat.quantity, UnknownQuantity)

    repeat = uncertain.description.sequence[0]
    assert repeat.is_quantity_uncertain is True
    assert repeat.is_quantity_known is False
    assert repeat.is_quantity_unknown is False
    assert isinstance(repeat.quantity, UncertainQuantity)
    assert repeat.quantity.lo == 19
    assert repeat.quantity.hi == 23

    repeat = omitted_unit.description.edit.sequence[0]
    assert repeat.unit is None
    assert repeat.is_unit_known is False
    assert repeat.is_quantity_known is True
    assert isinstance(repeat.quantity, KnownQuantity)


def test_parses_rna_special_outcomes_and_helper_views():
    produced = parse_hgvs("NM_004006.3:r.(4072_5145del)")
    no_change = parse_hgvs("NM_004006.3:r.=")
    predicted_no_change = parse_hgvs("NM_004006.3:r.(=)")
    not_produced = parse_hgvs("NM_004006.3:r.0")
    predicted_not_produced = parse_hgvs("NM_004006.3:r.0?")
    uncertain_splicing = parse_hgvs("NM_004006.3:r.spl")
    unknown = parse_hgvs("NM_004006.3:r.?")
    indeterminate = parse_hgvs("NM_004006.3:r.(?)")

    assert isinstance(produced.description, RnaProduced)
    assert produced.description.is_produced is True
    assert produced.description.is_no_change is False
    assert produced.description.is_not_produced is False
    assert produced.description.is_uncertain_splicing is False
    assert produced.description.is_unknown is False
    assert produced.description.is_indeterminate is False
    assert produced.description.is_predicted is True
    assert produced.description.certainty is OutcomeCertainty.PREDICTED
    assert isinstance(produced.description.edit, NucleotideDeletion)
    assert produced.description.edit.location.start.coordinate == 4072

    assert isinstance(no_change.description, RnaNoChange)
    assert no_change.description.is_no_change is True
    assert no_change.description.is_produced is False
    assert no_change.description.is_predicted is False
    assert no_change.description.certainty is OutcomeCertainty.CERTAIN

    assert isinstance(predicted_no_change.description, RnaNoChange)
    assert predicted_no_change.description.is_no_change is True
    assert predicted_no_change.description.is_predicted is True
    assert predicted_no_change.description.certainty is OutcomeCertainty.PREDICTED

    assert isinstance(not_produced.description, RnaNotProduced)
    assert not_produced.description.is_not_produced is True
    assert not_produced.description.is_no_change is False
    assert not_produced.description.is_predicted is False

    assert isinstance(predicted_not_produced.description, RnaNotProduced)
    assert predicted_not_produced.description.is_not_produced is True
    assert predicted_not_produced.description.is_predicted is True

    assert isinstance(uncertain_splicing.description, RnaUncertainSplicing)
    assert uncertain_splicing.description.is_uncertain_splicing is True
    assert uncertain_splicing.description.is_produced is False

    assert isinstance(unknown.description, RnaUnknown)
    assert unknown.description.is_unknown is True
    assert unknown.description.is_indeterminate is False

    assert isinstance(indeterminate.description, RnaIndeterminate)
    assert indeterminate.description.is_indeterminate is True
    assert indeterminate.description.is_unknown is False


def test_parses_nucleotide_allele_variants():
    cis = parse_hgvs("NC_000001.11:g.[123G>A;345del]")
    trans = parse_hgvs("NM_004006.3:r.[123c>a];[345del]")
    uncertain = parse_hgvs("NC_000001.11:g.123G>A(;)345del")
    unchanged = parse_hgvs("NM_004006.2:c.[2376G>C];[2376=]")
    mixed = parse_hgvs("NC_000001.11:g.[123G>A];[345del](;)789dup")

    assert isinstance(cis.description, AlleleVariant)
    assert len(cis.description.allele_one.variants) == 2
    assert cis.description.allele_two is None
    assert cis.description.phase is None
    assert cis.description.unphased == ()
    assert cis.description.allele_one.variants[0].reference == "G"
    assert cis.description.allele_one.variants[0].alternate == "A"
    assert (
        cis.description.allele_one.variants[1].location.start.coordinate == 345
    )
    assert isinstance(cis.description.allele_one.variants[1], NucleotideDeletion)

    assert isinstance(trans.description, AlleleVariant)
    assert len(trans.description.allele_one.variants) == 1
    assert trans.description.phase is AllelePhase.TRANS
    assert trans.description.allele_two is not None
    assert len(trans.description.allele_two.variants) == 1
    assert trans.description.unphased == ()
    assert (
        trans.description.allele_two.variants[0].edit.location.start.coordinate
        == 345
    )

    assert isinstance(uncertain.description, AlleleVariant)
    assert len(uncertain.description.allele_one.variants) == 1
    assert uncertain.description.phase is AllelePhase.UNCERTAIN
    assert uncertain.description.allele_two is not None
    assert (
        uncertain.description.allele_two.variants[0].location.start.coordinate
        == 345
    )
    assert uncertain.description.unphased == ()
    assert isinstance(uncertain.description.allele_two.variants[0], NucleotideDeletion)

    assert isinstance(unchanged.description, AlleleVariant)
    assert len(unchanged.description.allele_one.variants) == 1
    assert unchanged.description.phase is AllelePhase.TRANS
    assert unchanged.description.allele_two is not None
    assert (
        unchanged.description.allele_two.variants[0].location.start.coordinate
        == 2376
    )
    assert isinstance(unchanged.description.allele_two.variants[0], NucleotideNoChange)

    assert isinstance(mixed.description, AlleleVariant)
    assert len(mixed.description.allele_one.variants) == 1
    assert mixed.description.phase is AllelePhase.TRANS
    assert mixed.description.allele_two is not None
    assert len(mixed.description.unphased) == 1
    assert (
        mixed.description.allele_two.variants[0].location.start.coordinate
        == 345
    )
    assert mixed.description.unphased[0].location.start.coordinate == 789


def test_reports_nucleotide_allele_helper_views():
    cis = parse_hgvs("NC_000001.11:g.[123G>A;345del]")
    trans = parse_hgvs("NM_004006.3:r.[123c>a];[345del]")
    uncertain = parse_hgvs("NC_000001.11:g.123G>A(;)345del")
    mixed = parse_hgvs("NC_000001.11:g.[123G>A];[345del](;)789dup")

    assert cis.description.phased_alleles is None
    assert cis.description.unphased == ()
    assert len(cis.description.allele_one.variants) == 2
    assert cis.description.allele_two is None
    assert len(tuple(cis.description.allele_one)) == 2

    trans_pair = trans.description.phased_alleles
    assert trans_pair is not None
    assert len(trans_pair[0].variants) == 1
    assert len(trans_pair[1].variants) == 1
    assert trans.description.unphased == ()
    assert len(trans.description.allele_one.variants) == 1
    assert trans.description.allele_two is not None
    assert (
        trans.description.allele_two.variants[0].edit.location.start.coordinate
        == 345
    )

    assert uncertain.description.phased_alleles is None
    assert uncertain.description.unphased == ()
    assert len(uncertain.description.allele_one.variants) == 1
    assert uncertain.description.allele_two is not None
    assert (
        uncertain.description.allele_two.variants[0].location.start.coordinate
        == 345
    )

    mixed_pair = mixed.description.phased_alleles
    assert mixed_pair is not None
    assert len(mixed.description.unphased) == 1
    assert len(mixed.description.allele_one.variants) == 1
    assert mixed.description.allele_two is not None
    assert (
        mixed.description.allele_two.variants[0].location.start.coordinate
        == 345
    )
    assert mixed.description.unphased[0].location.start.coordinate == 789


def test_reports_derived_and_alternative_allele_forms():
    derived = parse_hgvs("NP_003997.1:p.[Ser68Arg,Asn594del]")
    alternative = parse_hgvs("NP_003997.1:p.[Ser68Arg]^[Asn594del]")
    regular = parse_hgvs("NP_003997.1:p.[Ser68Arg];[Ser68=]")
    mixed = parse_hgvs("NC_000001.11:g.[123G>A];[345del](;)789dup")

    assert isinstance(derived.description, DerivedAlleleForm)
    assert derived.description.is_derived is True
    assert derived.description.is_single is False
    assert derived.description.is_alternative is False
    assert len(derived.description.outcomes) == 2
    assert derived.description.outcomes[0].edit.to == "Arg"
    assert isinstance(derived.description.outcomes[1].edit, ProteinDeletion)

    assert isinstance(alternative.description, AlternativeAlleleForm)
    assert alternative.description.is_alternative is True
    assert alternative.description.is_single is False
    assert alternative.description.is_derived is False
    assert len(alternative.description.alternatives) == 2
    assert (
        alternative.description.alternatives[
            0
        ].allele_one.variants[0].edit.to
        == "Arg"
    )
    assert isinstance(
        alternative.description.alternatives[1].allele_one.variants[0].edit,
        ProteinDeletion,
    )

    assert isinstance(regular.description, AlleleVariant)
    assert regular.description.is_single is True
    assert regular.description.is_derived is False
    assert regular.description.is_alternative is False
    assert regular.description.has_second_allele is True
    assert regular.description.has_unphased_variants is False

    assert isinstance(mixed.description, AlleleVariant)
    assert mixed.description.has_second_allele is True
    assert mixed.description.has_unphased_variants is True


def test_rejects_malformed_nucleotide_allele_variants():
    cases = [
        "NC_000001.11:g.[123G>A](;)345del",
        "NC_000001.11:g.123G>A(;)[345del]",
        "NC_000001.11:g.[123G>A](;)[345del]",
        "NC_000001.11:g.[123G>A;;345del]",
        "NC_000001.11:g.[123G>A](;)",
        "NC_000001.11:g.[123G>A][345del]",
        "NC_000001.11:g.[123G>A;]",
        "NM_004006.3:r.;[123c>a]",
    ]

    for input_value in cases:
        with pytest.raises(TinyHGVSError) as exc_info:
            parse_hgvs(input_value)

        assert exc_info.value.code == "invalid.syntax"


def test_parses_protein_allele_variants():
    single = parse_hgvs("p.[Ser73Arg]")
    cis = parse_hgvs("NP_003997.1:p.[Ser68Arg;Asn594del]")
    trans = parse_hgvs("NP_003997.1:p.[Ser68Arg];[Ser68=]")
    uncertain = parse_hgvs("NP_003997.1:p.(Ser73Arg)(;)(Asn103del)")
    absent = parse_hgvs("p.[Ser86Arg];[0]")
    mixed = parse_hgvs("p.[Phe233Leu;(Cys690Trp)]")
    whole_predicted = parse_hgvs("NP_003997.1:p.[(Ser68Arg;Asn594del)]")
    predicted_with_no_change = parse_hgvs("NP_003997.1:p.[(Ser68Arg;Ser68=)]")
    range_no_change = parse_hgvs("p.[Ser68_Arg70dup];[Ser68_Arg70=]")

    assert isinstance(single.description, AlleleVariant)
    assert len(single.description.allele_one.variants) == 1
    assert single.description.allele_two is None
    assert single.description.phased_alleles is None

    assert isinstance(cis.description, AlleleVariant)
    assert len(cis.description.allele_one.variants) == 2
    assert cis.description.allele_two is None
    assert cis.description.phased_alleles is None
    assert cis.description.unphased == ()
    assert cis.description.allele_one.variants[0].is_predicted is False
    assert isinstance(cis.description.allele_one.variants[0], ProteinProduced)
    assert (
        cis.description.allele_one.variants[0].edit.location.start.residue
        == "Ser"
    )
    assert isinstance(cis.description.allele_one.variants[0].edit, ProteinSubstitution)
    assert cis.description.allele_one.variants[0].edit.to == "Arg"

    assert isinstance(trans.description, AlleleVariant)
    assert trans.description.phase is AllelePhase.TRANS
    assert trans.description.allele_two is not None
    assert trans.description.phased_alleles is not None
    assert isinstance(trans.description.allele_two.variants[0], ProteinNoChange)
    assert (
        trans.description.allele_two.variants[0].location.start.residue
        == "Ser"
    )
    assert trans.description.allele_two.variants[0].location.start.ordinal == 68
    assert trans.description.allele_two.variants[0].location.end is None

    assert isinstance(uncertain.description, AlleleVariant)
    assert uncertain.description.phase is AllelePhase.UNCERTAIN
    assert uncertain.description.allele_two is not None
    assert uncertain.description.phased_alleles is None
    assert uncertain.description.unphased == ()
    assert uncertain.description.allele_one.variants[0].is_predicted is True
    assert uncertain.description.allele_two.variants[0].is_predicted is True

    assert isinstance(absent.description, AlleleVariant)
    assert absent.description.phase is AllelePhase.TRANS
    assert absent.description.allele_two is not None
    assert isinstance(absent.description.allele_two.variants[0], ProteinNotProduced)

    assert isinstance(mixed.description, AlleleVariant)
    assert len(mixed.description.allele_one.variants) == 2
    assert mixed.description.allele_one.variants[0].is_predicted is False
    assert mixed.description.allele_one.variants[1].is_predicted is True

    assert isinstance(whole_predicted.description, AlleleVariant)
    assert len(whole_predicted.description.allele_one.variants) == 2
    assert all(
        variant.is_predicted
        for variant in whole_predicted.description.allele_one.variants
    )

    assert isinstance(predicted_with_no_change.description, AlleleVariant)
    assert len(predicted_with_no_change.description.allele_one.variants) == 2
    assert all(
        variant.is_predicted
        for variant in predicted_with_no_change.description.allele_one.variants
    )
    second_predicted = predicted_with_no_change.description.allele_one.variants[1]
    assert isinstance(second_predicted, ProteinNoChange)
    assert second_predicted.location.start.residue == "Ser"
    assert second_predicted.location.start.ordinal == 68
    assert second_predicted.certainty is OutcomeCertainty.PREDICTED

    assert isinstance(range_no_change.description, AlleleVariant)
    assert range_no_change.description.allele_two is not None
    second_range = range_no_change.description.allele_two.variants[0]
    assert isinstance(second_range, ProteinNoChange)
    assert second_range.location.start.residue == "Ser"
    assert second_range.location.start.ordinal == 68
    assert second_range.location.end is not None
    assert second_range.location.end.residue == "Arg"
    assert second_range.location.end.ordinal == 70


def test_reports_protein_allele_helper_views():
    single = parse_hgvs("p.[Ser73Arg]")
    trans = parse_hgvs("NP_003997.1:p.[Ser68Arg];[Ser68=]")
    uncertain = parse_hgvs("NP_003997.1:p.(Ser73Arg)(;)(Asn103del)")

    assert single.description.phased_alleles is None
    assert single.description.unphased == ()

    trans_pair = trans.description.phased_alleles
    assert trans_pair is not None
    assert len(trans_pair[0].variants) == 1
    assert len(trans_pair[1].variants) == 1
    assert trans.description.unphased == ()

    assert uncertain.description.phased_alleles is None
    assert uncertain.description.unphased == ()
    assert uncertain.description.allele_two is not None
    assert (
        uncertain.description.allele_two.variants[
            0
        ].edit.location.start.residue
        == "Asn"
    )


def test_rejects_malformed_protein_allele_variants():
    cases = [
        "p.([Ser68Arg;Asn594del])",
        "p.([Ser68Arg];[Ser68Arg])",
        "p.[=]",
        "p.[(=)]",
        "p.[Ser68Arg;=]",
        "p.[(Ser68Arg;=)]",
        "p.[Ser68Arg];[=]",
        "p.[Ser73Arg];[]",
        "p.[Ser68Arg](;)Asn594del",
        "p.[Ser73Arg+p.Asn103del]",
        "p.[Ser73Arg;p.Asn103del]",
    ]

    for input_value in cases:
        with pytest.raises(TinyHGVSError) as exc_info:
            parse_hgvs(input_value)

        assert exc_info.value.code == "invalid.syntax"


def test_parses_protein_substitution_and_no_change_variants():
    substitution = parse_hgvs("NP_003997.1:p.Trp24Ter")
    no_change = parse_hgvs("NP_003997.1:p.Cys188=")
    predicted_no_change = parse_hgvs("NP_003997.1:p.(Cys188=)")
    predicted_whole_no_change = parse_hgvs("NP_003997.1:p.(=)")

    assert isinstance(substitution.description, ProteinProduced)
    assert substitution.description.is_predicted is False
    assert isinstance(substitution.description.edit, ProteinSubstitution)
    assert substitution.description.edit.location.start.residue == "Trp"
    assert substitution.description.edit.location.start.ordinal == 24
    assert substitution.description.edit.to == "Ter"

    assert isinstance(no_change.description, ProteinNoChange)
    assert no_change.description.location.start.residue == "Cys"
    assert no_change.description.location.start.ordinal == 188
    assert no_change.description.location.end is None
    assert no_change.description.certainty is OutcomeCertainty.CERTAIN
    assert no_change.description.is_predicted is False

    assert isinstance(predicted_no_change.description, ProteinNoChange)
    assert predicted_no_change.description.location.start.residue == "Cys"
    assert predicted_no_change.description.location.start.ordinal == 188
    assert predicted_no_change.description.location.end is None
    assert predicted_no_change.description.certainty is OutcomeCertainty.PREDICTED
    assert predicted_no_change.description.is_predicted is True

    assert isinstance(predicted_whole_no_change.description, ProteinNoChange)
    assert predicted_whole_no_change.description.location is None
    assert (
        predicted_whole_no_change.description.certainty
        is OutcomeCertainty.PREDICTED
    )
    assert predicted_whole_no_change.description.is_predicted is True


def test_rejects_malformed_protein_no_change_variants():
    cases = [
        "NP_003997.1:p.Cys188(=)",
        "NP_003997.1:p.(Cys188(=))",
    ]

    for input_value in cases:
        with pytest.raises(TinyHGVSError) as exc_info:
            parse_hgvs(input_value)

        assert exc_info.value.code == "invalid.syntax"


def test_parses_uncertain_protein_locations():
    variant = parse_hgvs("NP_003997.1:p.(Ala123_Pro131)Ter")

    assert isinstance(variant.description, ProteinProduced)
    assert isinstance(variant.description.edit, ProteinSubstitution)
    location = variant.description.edit.location
    assert isinstance(location, UncertainLocation)
    assert location.is_interval is True
    assert location.is_position is False
    assert location.start.start.residue == "Ala"
    assert location.start.start.ordinal == 123
    assert location.start.end is not None
    assert location.start.end.residue == "Pro"
    assert location.start.end.ordinal == 131
    assert location.end is None
    assert variant.description.edit.to == "Ter"


def test_parses_protein_unknown_and_predicted_effects():
    unknown = parse_hgvs("NP_003997.1:p.?")
    predicted = parse_hgvs("NP_003997.1:p.(Trp24Ter)")
    absent = parse_hgvs("LRG_199p1:p.0")

    assert isinstance(unknown.description, ProteinUnknown)

    assert isinstance(predicted.description, ProteinProduced)
    assert predicted.description.is_predicted is True
    assert predicted.description.edit.location.start.residue == "Trp"
    assert predicted.description.edit.location.start.ordinal == 24

    assert isinstance(absent.description, ProteinNotProduced)


def test_parses_protein_deletion_duplication_insertion_and_delins_variants():
    deletion = parse_hgvs("NP_003997.2:p.Lys23_Val25del")
    duplication = parse_hgvs("NP_003997.2:p.Val7dup")
    insertion = parse_hgvs("p.Lys2_Gly3insGlnSerLys")
    delins = parse_hgvs("p.Cys28delinsTrpVal")

    assert isinstance(deletion.description.edit, ProteinDeletion)
    assert deletion.description.edit.location.start.residue == "Lys"
    assert deletion.description.edit.location.end is not None
    assert deletion.description.edit.location.end.residue == "Val"

    assert isinstance(duplication.description.edit, ProteinDuplication)

    assert isinstance(insertion.description.edit, KnownProteinInsertion)
    assert insertion.description.edit.sequence == (
        "Gln",
        "Ser",
        "Lys",
    )

    assert isinstance(delins.description.edit, ProteinDeletionInsertion)
    assert delins.description.edit.sequence == ("Trp", "Val")


def test_parses_protein_insertion_content_forms():
    known = parse_hgvs("NP_003997.1:p.Lys2_Gly3insGlnSerLys")
    unknown_single = parse_hgvs("NP_003997.1:p.Ser332_Ser333insXaa")
    unknown_count = parse_hgvs("NP_003997.1:p.Arg78_Gly79insXaa[23]")
    terminating = parse_hgvs("NP_003997.1:p.Gln746_Lys747ins*63")

    edit = known.description.edit
    assert isinstance(edit, KnownProteinInsertion)
    assert edit.is_insertion is True
    assert edit.is_deletion is False
    assert edit.location.start.residue == "Lys"
    assert edit.location.end is not None
    assert edit.location.end.residue == "Gly"
    assert edit.sequence == ("Gln", "Ser", "Lys")

    edit = unknown_single.description.edit
    assert isinstance(edit, UnknownProteinInsertion)
    assert edit.is_insertion is True
    assert edit.count == 1
    assert edit.location.start.residue == "Ser"
    assert edit.location.end is not None
    assert edit.location.end.residue == "Ser"

    edit = unknown_count.description.edit
    assert isinstance(edit, UnknownProteinInsertion)
    assert edit.is_insertion is True
    assert edit.count == 23

    edit = terminating.description.edit
    assert isinstance(edit, TerminatingProteinInsertion)
    assert edit.is_insertion is True
    assert edit.ordinal == 63
    assert edit.location.start.residue == "Gln"
    assert edit.location.end is not None
    assert edit.location.end.residue == "Lys"


def test_parses_protein_uncertain_consequence_forms():
    substitution = parse_hgvs("NP_003997.1:p.(Gly719Ala^Ser)")
    frameshift = parse_hgvs(
        "NP_003997.1:p.(Gly23GlufsTer7^Gly23CysfsTer26)"
    )

    assert isinstance(substitution.description, ProteinProduced)
    assert substitution.description.is_produced is True
    assert substitution.description.has_alternatives is False
    assert substitution.description.is_predicted is True
    assert substitution.description.certainty is OutcomeCertainty.PREDICTED
    assert isinstance(substitution.description.edit, ProteinSubstitution)
    assert substitution.description.edit.has_alternatives is True
    assert substitution.description.edit.to == ("Ala", "Ser")

    assert isinstance(frameshift.description, ProteinProducedAlternatives)
    assert frameshift.description.is_produced is True
    assert frameshift.description.is_no_change is False
    assert frameshift.description.is_unknown is False
    assert frameshift.description.has_alternatives is True
    assert frameshift.description.is_predicted is True
    assert frameshift.description.certainty is OutcomeCertainty.PREDICTED
    assert len(frameshift.description.edits) == 2
    assert all(edit.is_frameshift for edit in frameshift.description.edits)
    assert frameshift.description.edits[0].location.start.residue == "Gly"
    assert frameshift.description.edits[0].to_residue == "Glu"
    assert frameshift.description.edits[1].to_residue == "Cys"


def test_parses_protein_repeat_variants():
    repeat = parse_hgvs("NP_0123456.1:p.Arg65_Ser67[12]")

    assert isinstance(repeat.description.edit, ProteinRepeat)
    assert repeat.description.edit.location.start.residue == "Arg"
    assert repeat.description.edit.location.start.ordinal == 65
    assert repeat.description.edit.location.end is not None
    assert repeat.description.edit.location.end.residue == "Ser"
    assert repeat.description.edit.location.end.ordinal == 67
    assert repeat.description.edit.repeat.quantity.count == 12


def test_parses_protein_frameshift_variants():
    short = parse_hgvs("NP_0123456.1:p.Arg97fs")
    long = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
    symbolic_stop = parse_hgvs("NP_0123456.1:p.Arg97Profs*23")
    unknown_stop = parse_hgvs("NP_0123456.1:p.Ile327Argfs*?")
    unknown_stop_ter = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer?")
    predicted = parse_hgvs("p.(Arg97fs)")

    assert isinstance(short.description.edit, ProteinFrameshift)
    assert short.description.edit.location.start.residue == "Arg"
    assert short.description.edit.location.start.ordinal == 97
    assert short.description.edit.to_residue is None
    assert isinstance(short.description.edit.stop, OmittedProteinFrameshiftStop)

    assert isinstance(long.description.edit, ProteinFrameshift)
    assert long.description.edit.to_residue == "Pro"
    assert isinstance(long.description.edit.stop, KnownProteinFrameshiftStop)
    assert long.description.edit.stop.ordinal == 23

    assert isinstance(symbolic_stop.description.edit, ProteinFrameshift)
    assert symbolic_stop.description.edit.to_residue == "Pro"
    assert isinstance(symbolic_stop.description.edit.stop, KnownProteinFrameshiftStop)
    assert symbolic_stop.description.edit.stop.ordinal == 23

    assert isinstance(unknown_stop.description.edit, ProteinFrameshift)
    assert unknown_stop.description.edit.location.start.residue == "Ile"
    assert unknown_stop.description.edit.location.start.ordinal == 327
    assert unknown_stop.description.edit.to_residue == "Arg"
    assert isinstance(unknown_stop.description.edit.stop, UnknownProteinFrameshiftStop)

    assert isinstance(unknown_stop_ter.description.edit, ProteinFrameshift)
    assert unknown_stop_ter.description.edit.to_residue == "Pro"
    assert isinstance(
        unknown_stop_ter.description.edit.stop, UnknownProteinFrameshiftStop
    )

    assert isinstance(predicted.description.edit, ProteinFrameshift)
    assert predicted.description.is_predicted is True
    assert predicted.description.edit.to_residue is None
    assert isinstance(predicted.description.edit.stop, OmittedProteinFrameshiftStop)


def test_reports_protein_frameshift_stop_helper_views():
    omitted = parse_hgvs("NP_0123456.1:p.Arg97fs")
    known = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
    unknown = parse_hgvs("NP_0123456.1:p.Ile327Argfs*?")

    stop = omitted.description.edit.stop
    assert isinstance(stop, OmittedProteinFrameshiftStop)
    assert stop.is_omitted is True
    assert stop.is_known is False
    assert stop.is_unknown is False

    stop = known.description.edit.stop
    assert isinstance(stop, KnownProteinFrameshiftStop)
    assert stop.is_known is True
    assert stop.is_omitted is False
    assert stop.is_unknown is False
    assert stop.ordinal == 23

    stop = unknown.description.edit.stop
    assert isinstance(stop, UnknownProteinFrameshiftStop)
    assert stop.is_unknown is True
    assert stop.is_omitted is False
    assert stop.is_known is False


def test_reports_protein_edit_and_outcome_helper_views():
    substitution = parse_hgvs("NP_003997.1:p.Trp24Ter")
    deletion = parse_hgvs("NP_003997.2:p.Lys23_Val25del")
    duplication = parse_hgvs("NP_003997.2:p.Val7dup")
    repeat = parse_hgvs("NP_0123456.1:p.Arg65_Ser67[12]")
    extension = parse_hgvs("NP_003997.2:p.Met1ext-5")
    frameshift = parse_hgvs("NP_0123456.1:p.Arg97fs")
    insertion = parse_hgvs("NP_003997.1:p.Lys2_Gly3insGlnSerLys")
    delins = parse_hgvs("NP_003997.1:p.Cys28delinsTrpVal")
    no_change = parse_hgvs("NP_003997.1:p.(Cys188=)")
    not_produced = parse_hgvs("LRG_199p1:p.0")
    unknown = parse_hgvs("NP_003997.1:p.?")

    edit = substitution.description.edit
    assert edit.is_substitution is True
    assert edit.is_deletion is False
    assert edit.is_duplication is False
    assert edit.is_repeat is False
    assert edit.is_extension is False
    assert edit.is_frameshift is False
    assert edit.is_insertion is False
    assert edit.is_deletion_insertion is False
    assert edit.has_alternatives is False

    assert deletion.description.edit.is_deletion is True
    assert deletion.description.edit.is_substitution is False
    assert duplication.description.edit.is_duplication is True
    assert repeat.description.edit.is_repeat is True
    assert extension.description.edit.is_extension is True
    assert frameshift.description.edit.is_frameshift is True
    assert frameshift.description.edit.has_alternative_residues is False
    assert insertion.description.edit.is_insertion is True
    assert delins.description.edit.is_deletion_insertion is True

    assert substitution.description.is_produced is True
    assert substitution.description.is_no_change is False
    assert substitution.description.is_not_produced is False
    assert substitution.description.is_unknown is False

    assert isinstance(no_change.description, ProteinNoChange)
    assert no_change.description.is_no_change is True
    assert no_change.description.is_predicted is True
    assert no_change.description.certainty is OutcomeCertainty.PREDICTED

    assert isinstance(not_produced.description, ProteinNotProduced)
    assert not_produced.description.is_not_produced is True
    assert not_produced.description.is_produced is False
    assert not_produced.description.is_predicted is False

    assert isinstance(unknown.description, ProteinUnknown)
    assert unknown.description.is_unknown is True
    assert unknown.description.is_produced is False


def test_parses_protein_extension_variants():
    n_terminal = parse_hgvs("NP_003997.2:p.Met1ext-5")
    predicted_n_terminal = parse_hgvs("p.(Met1ext-8)")
    c_terminal = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")
    c_terminal_symbolic = parse_hgvs("p.*110Glnext*17")
    unknown_stop = parse_hgvs("p.Ter327ArgextTer?")
    unknown_stop_symbolic = parse_hgvs("p.*327Argext*?")

    assert isinstance(n_terminal.description.edit, ProteinExtension)
    assert n_terminal.description.edit.location.start.residue == "Met"
    assert n_terminal.description.edit.location.start.ordinal == 1
    assert n_terminal.description.edit.to_terminal is ProteinExtensionTerminal.N
    assert n_terminal.description.edit.to_residue is None
    assert n_terminal.description.edit.terminal_ordinal == -5

    assert isinstance(predicted_n_terminal.description.edit, ProteinExtension)
    assert predicted_n_terminal.description.is_predicted is True
    assert (
        predicted_n_terminal.description.edit.to_terminal
        is ProteinExtensionTerminal.N
    )
    assert predicted_n_terminal.description.edit.terminal_ordinal == -8

    assert isinstance(c_terminal.description.edit, ProteinExtension)
    assert c_terminal.description.edit.location.start.residue == "Ter"
    assert c_terminal.description.edit.location.start.ordinal == 110
    assert c_terminal.description.edit.to_terminal is ProteinExtensionTerminal.C
    assert c_terminal.description.edit.to_residue == "Gln"
    assert c_terminal.description.edit.terminal_ordinal == 17

    assert isinstance(c_terminal_symbolic.description.edit, ProteinExtension)
    assert (
        c_terminal_symbolic.description.edit.location.start.residue == "Ter"
    )
    assert c_terminal_symbolic.description.edit.location.start.ordinal == 110
    assert (
        c_terminal_symbolic.description.edit.to_terminal
        is ProteinExtensionTerminal.C
    )
    assert c_terminal_symbolic.description.edit.to_residue == "Gln"
    assert c_terminal_symbolic.description.edit.terminal_ordinal == 17

    assert isinstance(unknown_stop.description.edit, ProteinExtension)
    assert unknown_stop.description.edit.location.start.residue == "Ter"
    assert unknown_stop.description.edit.location.start.ordinal == 327
    assert (
        unknown_stop.description.edit.to_terminal
        is ProteinExtensionTerminal.C
    )
    assert unknown_stop.description.edit.to_residue == "Arg"
    assert unknown_stop.description.edit.terminal_ordinal is None

    assert isinstance(unknown_stop_symbolic.description.edit, ProteinExtension)
    assert (
        unknown_stop_symbolic.description.edit.location.start.residue
        == "Ter"
    )
    assert (
        unknown_stop_symbolic.description.edit.location.start.ordinal == 327
    )
    assert (
        unknown_stop_symbolic.description.edit.to_terminal
        is ProteinExtensionTerminal.C
    )
    assert unknown_stop_symbolic.description.edit.to_residue == "Arg"
    assert unknown_stop_symbolic.description.edit.terminal_ordinal is None


def test_reports_intronic_and_utr_coordinate_properties_from_parsed_variants():
    intronic = parse_hgvs("NM_004006.2:c.93+1G>T").description.location.start
    five_prime_intronic = parse_hgvs(
        "NM_001385026.1:c.-106+2T>A"
    ).description.location.start
    five_prime_utr = parse_hgvs(
        "NM_007373.4:c.-81C>T"
    ).description.location.start
    three_prime_intronic = parse_hgvs(
        "NM_001272071.2:c.*639-1G>A"
    ).description.location.start
    three_prime_utr = parse_hgvs(
        "NM_001272071.2:c.*1C>T"
    ).description.location.start

    assert intronic.is_intronic is True
    assert intronic.is_cds_start_anchored is False
    assert intronic.is_cds_end_anchored is False
    assert intronic.is_five_prime_utr is False
    assert intronic.is_three_prime_utr is False

    assert five_prime_intronic.is_intronic is True
    assert five_prime_intronic.is_cds_start_anchored is True
    assert five_prime_intronic.is_cds_end_anchored is False
    assert five_prime_intronic.is_five_prime_utr is False
    assert five_prime_intronic.is_three_prime_utr is False

    assert five_prime_utr.is_intronic is False
    assert five_prime_utr.is_cds_start_anchored is True
    assert five_prime_utr.is_cds_end_anchored is False
    assert five_prime_utr.is_five_prime_utr is True
    assert five_prime_utr.is_three_prime_utr is False

    assert three_prime_intronic.is_intronic is True
    assert three_prime_intronic.is_cds_start_anchored is False
    assert three_prime_intronic.is_cds_end_anchored is True
    assert three_prime_intronic.is_five_prime_utr is False
    assert three_prime_intronic.is_three_prime_utr is False

    assert three_prime_utr.is_intronic is False
    assert three_prime_utr.is_cds_start_anchored is False
    assert three_prime_utr.is_cds_end_anchored is True
    assert three_prime_utr.is_five_prime_utr is False
    assert three_prime_utr.is_three_prime_utr is True


def test_rejects_malformed_nucleotide_coordinate_model_states():
    known = NucleotideCoordinate(NucleotideAnchor.ABSOLUTE, 123, 0)
    unknown = NucleotideCoordinate(None, None, None)

    assert known.is_known is True
    assert known.is_unknown is False
    assert unknown.is_unknown is True
    assert unknown.is_known is False

    with pytest.raises(ValueError):
        NucleotideCoordinate(NucleotideAnchor.ABSOLUTE, 123, None)

    with pytest.raises(ValueError):
        NucleotideCoordinate(None, 123, 0)


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
        "NM_001385026.1:c.-106+T>A",
        "NM_001385026.1:c.-106++2T>A",
        "NM_001272071.2:c.*639--1G>A",
        "NM_001272071.2:c.*24-12888_+5del",
        "NM_001385026.1:c.-0+2A>G",
        "NM_001272071.2:c.*0-1G>A",
    ],
)
def test_rejects_malformed_cdna_offset_anchor_variants(example: str):
    with pytest.raises(TinyHGVSError) as exc_info:
        parse_hgvs(example)

    assert exc_info.value.code == "invalid.syntax"
    assert exc_info.value.kind is ParseHgvsErrorKind.INVALID_SYNTAX


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


@pytest.mark.parametrize(
    ("example", "code", "kind", "fragment"),
    [
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
            "NM_002354.2:r.-358_555::NM_000251.2:r.212_*279",
            "unsupported.rna_adjoined_transcript",
            ParseHgvsErrorKind.UNSUPPORTED_SYNTAX,
            "::",
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
