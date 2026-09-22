"""Nucleotide-focused Python model types for parsed HGVS variants.

Type Aliases:
    NucleotideSequenceItem: Tagged union for supported inserted or replacement
        nucleotide sequence models.
    NucleotideEdit: Tagged union for supported nucleotide edit models.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Literal, TypeAlias

from .shared import (
    Allele,
    AllelePhase,
    AlleleVariant,
    CoordinateSystem,
    Interval,
    KnownLocation,
    Location,
    ReferenceSpec,
    Repeat,
)


class NucleotideAnchor(str, Enum):
    """HGVS reference point used to interpret a nucleotide position.

    Attributes:
        ABSOLUTE: Coordinate is read directly on the named reference sequence.
        RELATIVE_CDS_START: Coordinate is read relative to the CDS start site.
        RELATIVE_CDS_END: Coordinate is read relative to the CDS end site.

    Examples:
        An intronic splice-site substitution uses direct coordinates:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> variant.description.location.start.anchor
        <NucleotideAnchor.ABSOLUTE: 'absolute'>

        A 5' UTR substitution is anchored to the CDS start:
        >>> variant = parse_hgvs("NM_007373.4:c.-1C>T")
        >>> variant.description.location.start.anchor
        <NucleotideAnchor.RELATIVE_CDS_START: 'relative_cds_start'>

        A 3' UTR substitution is anchored to the CDS end:
        >>> variant = parse_hgvs("NM_001272071.2:c.*1C>T")
        >>> variant.description.location.start.anchor
        <NucleotideAnchor.RELATIVE_CDS_END: 'relative_cds_end'>
    """

    ABSOLUTE = "absolute"
    RELATIVE_CDS_START = "relative_cds_start"
    RELATIVE_CDS_END = "relative_cds_end"


class NucleotideCoordinateKind(str, Enum):
    """Known/unknown state for a nucleotide coordinate.

    Attributes:
        KNOWN: Coordinate has anchor, coordinate, and offset values.
        UNKNOWN: Coordinate is written as ``?``.
    """

    KNOWN = "known"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class NucleotideCoordinate:
    """A nucleotide coordinate or the explicit HGVS unknown coordinate `?`."""

    anchor: NucleotideAnchor | None
    coordinate: int | None
    offset: int | None

    @property
    def is_known(self) -> bool:
        return self.coordinate is not None

    @property
    def is_unknown(self) -> bool:
        return self.coordinate is None

    @property
    def is_intronic(self) -> bool:
        return self.offset not in (None, 0)

    @property
    def is_five_prime_utr(self) -> bool:
        return (
            self.anchor is NucleotideAnchor.RELATIVE_CDS_START
            and self.coordinate is not None
            and self.coordinate < 0
        )

    @property
    def is_three_prime_utr(self) -> bool:
        return self.anchor is NucleotideAnchor.RELATIVE_CDS_END


@dataclass(frozen=True, slots=True)
class NucleotdieNoChange:
    location: Location[NucleotideCoordinate]


@dataclass(frozen=True, slots=True)
class NucleotideSubstitution:
    location: Location[NucleotideCoordinate]
    reference: str
    alternate: str


@dataclass(frozen=True, slots=True)
class NucleotideDeletion:
    location: Location[NucleotideCoordinate]


@dataclass(frozen=True, slots=True)
class NucleotideDuplication:
    location: Location[NucleotideCoordinate]


@dataclass(frozen=True, slots=True)
class NucleotideInversion:
    location: Location[NucleotideCoordinate]


@dataclass(frozen=True, slots=True)
class NucleotideRepeat:
    location: Location[NucleotideCoordinate]
    sequence: tuple[Repeat, ...]


@dataclass(frozen=True, slots=True)
class NucleotideInsertion:
    location: Location[NucleotideCoordinate]
    sequence: tuple[NucleotideSequenceItem, ...]


@dataclass(frozen=True, slots=True)
class NucleotideDeletionInsertion:
    location: Location[NucleotideCoordinate]
    sequence: tuple[NucleotideSequenceItem, ...]


NucleotideEdit: TypeAlias = (
    NucleotdieNoChange
    | NucleotideSubstitution
    | NucleotideDeletion
    | NucleotideDuplication
    | NucleotideRepeat
    | NucleotideInsertion
    | NucleotideInversion
    | NucleotideDeletionInsertion
)


@dataclass(frozen=True, slots=True)
class LiteralSequence:
    value: str


@dataclass(frozen=True, slots=True)
class CopiedSequence:
    """Copied nucleotide sequence used in an insertion or deletion-insertion."""

    reference: ReferenceSpec | None
    coordinate_system: CoordinateSystem | None
    location: KnownLocation[NucleotideCoordinate]
    is_inverted: bool

    @property
    def is_from_same_reference(self) -> bool:
        return self.reference is None and self.coordinate_system is None


NucleotideSequenceItem: TypeAlias = LiteralSequence | Repeat | CopiedSequence

# @dataclass(frozen=True, slots=True)
# class CopiedSequenceItem:
#     """Copied nucleotide sequence used in an insertion or deletion-insertion.
#
#     Attributes:
#         source_reference: Source reference when the copied sequence comes from
#             a different accession. ``None`` means the same outer reference.
#         source_coordinate_system: Source coordinate system when it differs from
#             the outer variant. ``None`` means the same outer coordinate system.
#         source_location: Inclusive interval on the source reference.
#         is_inverted: Whether the copied sequence is inserted in reverse
#             orientation.
#
#     Examples:
#         A stretch of sequence from the same transcript is inserted in reverse
#         orientation:
#         >>> from tinyhgvs import parse_hgvs
#         >>> variant = parse_hgvs("NM_004006.2:c.849_850ins850_900inv")
#         >>> item = variant.description.edit.items[0]
#         >>> item.is_from_same_reference
#         True
#         >>> item.source_location.start.coordinate
#         850
#         >>> item.source_location.end.coordinate
#         900
#         >>> item.is_inverted
#         True
#
#         A copied sequence can also come from another chromosome:
#         >>> variant = parse_hgvs("NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]")
#         >>> item = variant.description.edit.items[0]
#         >>> item.source_reference.primary.id
#         'NC_000022.10'
#         >>> item.source_coordinate_system
#         <CoordinateSystem.GENOMIC: 'g'>
#     """
#
#     source_reference: ReferenceSpec | None
#     source_coordinate_system: CoordinateSystem | None
#     source_location: Interval[NucleotideCoordinate]
#     is_inverted: bool
#
#     @property
#     def is_from_same_reference(self) -> bool:
#         return (
#             self.source_reference is None
#             and self.source_coordinate_system is None
#         )


# @dataclass(frozen=True, slots=True)
# class LiteralSequenceItem:
#     """Model describing literal-base-type sequence edit component.
#
#     Attributes:
#         value: Nucleotide bases.
#
#     Examples:
#         A literal insertion of three nucleotides:
#         >>> from tinyhgvs import LiteralSequenceItem, parse_hgvs
#         >>> variant = parse_hgvs("NC_000023.10:g.32862923_32862924insCCT")
#         >>> item = variant.description.edit.items[0]
#         >>> isinstance(item, LiteralSequenceItem)
#         True
#         >>> item.value
#         'CCT'
#     """
#
#     value: str


@dataclass(frozen=True, slots=True)
class RepeatSequenceItem:
    """Model describing repeat-type sequence edit component.

    Attributes:
        unit: Repeat unit of nucleotide bases.
        count: Number of units being repeated.

    Examples:
        The insertion contains 100 copies of ``N``:
        >>> from tinyhgvs import RepeatSequenceItem, parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.32717298_32717299insN[100]")
        >>> item = variant.description.edit.items[0]
        >>> isinstance(item, RepeatSequenceItem)
        True
        >>> item.unit
        'N'
        >>> item.count
        100
    """

    unit: str
    count: int


@dataclass(frozen=True, slots=True)
class NucleotideRepeatBlock:
    """One repeat block/unit in a nucleotide repeat description.

    Attributes:
        count: Number of repeated units.
        unit: Literal base(s) repeat unit. None when repeat unit is described
            in the form of `location[count]`.
        location: Location per repeat block. None when repeat unit is described
            in the form of `unit[count]`.

    Examples:
        A literal 3bp bases repeat with 23 units:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000014.8:g.123CAG[23]")
        >>> block = variant.description.edit.blocks[0]
        >>> block.unit
        'CAG'
        >>> block.count
        23
        >>> block.location is None
        True

        A RNA repeat variant composed of consecutive repeat units, each
        described in the form `location[count]`, rather than `unit[count]`:
        a repetitive unit from a location:
        >>> variant = parse_hgvs("NM_004006.3:r.456_465[4]466_489[9]490_499[3]")
        >>> block = variant.description.edit.blocks[1]
        >>> block.unit is None
        True
        >>> block.location.start.coordinate
        466
        >>> block.location.end.coordinate
        489
    """

    count: int
    unit: str | None = None
    location: Interval[NucleotideCoordinate] | None = None


# NucleotideSequenceItem: TypeAlias = (
#     LiteralSequenceItem | RepeatSequenceItem | CopiedSequenceItem
# )
# """Tagged union for supported inserted or replacement nucleotide components:
#
# - [`LiteralSequenceItem`][tinyhgvs.models.nucleotide.LiteralSequenceItem]
# - [`RepeatSequenceItem`][tinyhgvs.models.nucleotide.RepeatSequenceItem]
# - [`CopiedSequenceItem`][tinyhgvs.models.nucleotide.CopiedSequenceItem]
# """


class NucleotideSequenceOmittedEdit(str, Enum):
    """Nucleotide edits whose altered sequence is not written explicitly.

    Attributes:
        NO_CHANGE: No nucleotide change, written as ``=``.
        DELETION: Deletion of the reference interval, written as ``del``.
        DUPLICATION: Duplication of the reference interval, written as ``dup``.
        INVERSION: Inversion of the reference interval, written as ``inv``.

    Examples:
        A coding DNA deletion:
        >>> from tinyhgvs import NucleotideSequenceOmittedEdit, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.5697del")
        >>> variant.description.edit is NucleotideSequenceOmittedEdit.DELETION
        True

        A genomic duplication:
        >>> variant = parse_hgvs("NC_000001.11:g.1234_2345dup")
        >>> variant.description.edit
        <NucleotideSequenceOmittedEdit.DUPLICATION: 'duplication'>
    """

    NO_CHANGE = "no_change"
    DELETION = "deletion"
    DUPLICATION = "duplication"
    INVERSION = "inversion"


@dataclass(frozen=True, slots=True)
class NucleotideSubstitutionEdit:
    """Model describing nucleotide substitution.

    Examples:
        A reference base ``C`` is substituted by ``A`` at the described
        location.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.33038255C>A")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.reference
        'C'
        >>> variant_edit.alternate
        'A'
        >>> variant_edit.kind
        'substitution'
    """

    reference: str
    alternate: str
    kind: Literal["substitution"] = field(init=False, default="substitution")


@dataclass(frozen=True, slots=True)
class NucleotideInsertionEdit:
    """Model describing nucleotide insertion.

    Attributes:
        items: Inserted sequence items in the order they appear in the HGVS
            expression.
        kind: Edit kind.

    Examples:
        Literal nucleotide insertion:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.32862923_32862924insCCT")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.items[0].value
        'CCT'

        A composite insertion can mix literal and copied sequence:
        >>> variant = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")
        >>> variant_edit = variant.description.edit
        >>> len(variant_edit.items)
        3
        >>> variant_edit.items[0].value
        'T'
        >>> variant_edit.items[1].is_from_same_reference
        True
        >>> variant_edit.items[2].value
        'AGGG'

        Insertion of unspecified repeated bases:
        >>> variant = parse_hgvs("NC_000023.10:g.32717298_32717299insN[100]")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.items[0].unit
        'N'
        >>> variant_edit.items[0].count
        100
    """

    items: tuple[NucleotideSequenceItem, ...]
    kind: Literal["insertion"] = field(init=False, default="insertion")


@dataclass(frozen=True, slots=True)
class NucleotideDeletionInsertionEdit:
    """Model describing nucleotide deletion-insertion.

    Attributes:
        items: Replacement sequence items in the order they appear in the HGVS
            expression.
        kind: Edit kind.

    Examples:
        A deleted interval is replaced by one literal sequence component.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("LRG_199t1:c.850_901delinsTTCCTCGATGCCTG")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.items[0].value
        'TTCCTCGATGCCTG'

        A deleted interval can be replaced by copied sequence from the same
        reference:
        >>> variant = parse_hgvs("NC_000022.10:g.42522624_42522669delins42536337_42536382")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.items[0].source_location.start.coordinate
        42536337

        A deleted interval is replaced by repeated unspecified bases.
        >>> variant = parse_hgvs("NM_004006.2:c.812_829delinsN[12]")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.items[0].unit
        'N'
        >>> variant_edit.items[0].count
        12
    """

    items: tuple[NucleotideSequenceItem, ...]
    kind: Literal["deletion_insertion"] = field(
        init=False,
        default="deletion_insertion",
    )


@dataclass(frozen=True, slots=True)
class NucleotideRepeatEdit:
    """Model describing a top-level nucleotide repeat variant.

    Attributes:
        blocks: Repeat blocks/units written in the HGVS description.
        kind: Edit kind.

    Examples:
        A DNA repeat variant with explicit repeat unit:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000014.8:g.123CAG[23]")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.blocks[0].unit
        'CAG'
        >>> variant_edit.blocks[0].count
        23

        A RNA repeat variant composed of consecutive blocks/units, each
        represented as a location span:
        >>> variant = parse_hgvs("NM_004006.3:r.456_465[4]466_489[9]490_499[3]")
        >>> variant_edit = variant.description.edit
        >>> len(variant_edit.blocks)
        3
        >>> variant_edit.blocks[2].count
        3
    """

    blocks: tuple[NucleotideRepeatBlock, ...]
    kind: Literal["repeat"] = field(init=False, default="repeat")


# NucleotideEdit: TypeAlias = (
#     NucleotideSequenceOmittedEdit
#     | NucleotideSubstitutionEdit
#     | NucleotideInsertionEdit
#     | NucleotideDeletionInsertionEdit
#     | NucleotideRepeatEdit
# )
# """Tagged union for supported nucleotide edit models:
#
# - [`NucleotideSequenceOmittedEdit`][tinyhgvs.models.nucleotide.NucleotideSequenceOmittedEdit]
# - [`NucleotideSubstitutionEdit`][tinyhgvs.models.nucleotide.NucleotideSubstitutionEdit]
# - [`NucleotideInsertionEdit`][tinyhgvs.models.nucleotide.NucleotideInsertionEdit]
# - [`NucleotideDeletionInsertionEdit`][tinyhgvs.models.nucleotide.NucleotideDeletionInsertionEdit]
# - [`NucleotideRepeatEdit`][tinyhgvs.models.nucleotide.NucleotideRepeatEdit]
# """


@dataclass(frozen=True, slots=True)
class NucleotideVariant:
    """Model describing a nucleotide-level variant.

    Attributes:
        location: [`Location`][tinyhgvs.models.shared.Location] where the
            nucleotide edit occurs.
        edit: Nucleotide edit applied at the location.

    Examples:
        A splice-site substitution is represented by a nucleotide location and
        a nucleotide substitution edit.
        >>> from tinyhgvs import NucleotideSubstitutionEdit, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> isinstance(variant.description.edit, NucleotideSubstitutionEdit)
        True
        >>> variant_description = variant.description
        >>> variant_description.location.start.coordinate
        357
        >>> variant_description.location.start.offset
        1
    """

    location: Location[NucleotideCoordinate]
    edit: NucleotideEdit


__all__ = [
    "Allele",
    "AllelePhase",
    "AlleleVariant",
    "CopiedSequenceItem",
    "NucleotideDeletionInsertionEdit",
    "NucleotideAnchor",
    "NucleotideCoordinate",
    "NucleotideCoordinateKind",
    "NucleotideEdit",
    "NucleotideInsertionEdit",
    "NucleotideRepeatBlock",
    "NucleotideRepeatEdit",
    "NucleotideSequenceItem",
    "NucleotideSequenceOmittedEdit",
    "NucleotideSubstitutionEdit",
    "NucleotideVariant",
    "LiteralSequenceItem",
    "RepeatSequenceItem",
]
