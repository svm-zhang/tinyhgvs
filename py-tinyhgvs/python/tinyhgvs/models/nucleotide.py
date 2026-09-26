"""Nucleotide-focused Python model types for parsed HGVS variants.

Type Aliases:
    NucleotideSequenceItem: Tagged union for supported inserted or replacement
        nucleotide sequence models.
    NucleotideEdit: Tagged union for supported nucleotide edit models.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TypeAlias

from .repeat import Repeat
from .shared import (
    CoordinateSystem,
    KnownLocation,
    Location,
    ReferenceSpec,
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


@dataclass(frozen=True, slots=True)
class NucleotideCoordinate:
    """A nucleotide coordinate or the explicit HGVS unknown coordinate `?`."""

    anchor: NucleotideAnchor | None
    coordinate: int | None
    offset: int | None

    def __post_init__(self) -> None:
        fields = (self.anchor, self.coordinate, self.offset)

        if all(field is None for field in fields):
            return

        if all(field is not None for field in fields):
            return

        # Reject malformed coordinate states that cannot represent valid HGVS
        # location pieces. Supported states are fully known coordinates like
        # 123 or 456 in (123_456), or fully unknown coordinates like ? in
        # (?_?) and (123_?). Mixed internal states with only some of anchor,
        # coordinate, and offset present cannot represent one valid boundary.
        raise ValueError(
            "NucleotideCoordinate must be either fully known or fully unknown"
        )

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
    def is_cds_start_anchored(self) -> bool:
        return self.anchor is NucleotideAnchor.RELATIVE_CDS_START

    @property
    def is_cds_end_anchored(self) -> bool:
        return self.anchor is NucleotideAnchor.RELATIVE_CDS_END

    @property
    def is_five_prime_utr(self) -> bool:
        return (
            self.is_cds_start_anchored
            and self.coordinate is not None
            and self.coordinate < 0
        )

    @property
    def is_three_prime_utr(self) -> bool:
        return self.is_cds_end_anchored


@dataclass(frozen=True, slots=True)
class NucleotideNoChange:
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
    NucleotideNoChange
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


__all__ = [
    "CopiedSequence",
    "NucleotideDeletionInsertion",
    "NucleotideAnchor",
    "NucleotideCoordinate",
    "NucleotideEdit",
    "NucleotideInsertion",
    "NucleotideDeletion",
    "NucleotideDuplication",
    "NucleotideRepeat",
    "NucleotideSequenceItem",
    "NucleotideSubstitution",
    "NucleotideNoChange",
    "LiteralSequence",
]
