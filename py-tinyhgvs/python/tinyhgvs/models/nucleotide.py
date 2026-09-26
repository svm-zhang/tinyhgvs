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
    """A nucleotide coordinate or the explicit HGVS unknown coordinate ``?``.

    Attributes:
        anchor: Reference point used to interpret the coordinate. ``None`` for
            unknown coordinates.
        coordinate: Primary HGVS coordinate as written. ``None`` for unknown
            coordinates.
        offset: Signed secondary displacement from the primary coordinate.
            ``None`` for unknown coordinates. Positive values move downstream
            and negative values move upstream.

    Examples:
        Duplication crossing an exon/intron border:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.11(NM_004006.2):c.260_264+48dup")
        >>> variant.description.location.start.coordinate
        260
        >>> variant.description.location.end.coordinate
        264
        >>> variant.description.location.end.offset
        48

        Upstream intronic substitution:
        >>> variant = parse_hgvs("NG_012232.1(NM_004006.2):c.264-2A>G")
        >>> variant.description.location.start.coordinate
        264
        >>> variant.description.location.start.offset
        -2

        5' UTR and 3' UTR coordinates keep their signed coordinate values:
        >>> five_prime = parse_hgvs("NM_007373.4:c.-1C>T")
        >>> five_prime.description.location.start.anchor
        <NucleotideAnchor.RELATIVE_CDS_START: 'relative_cds_start'>
        >>> five_prime.description.location.start.coordinate
        -1
        >>> three_prime = parse_hgvs("NM_001272071.2:c.*1C>T")
        >>> three_prime.description.location.start.anchor
        <NucleotideAnchor.RELATIVE_CDS_END: 'relative_cds_end'>
        >>> three_prime.description.location.start.coordinate
        1

        Unknown coordinates are used by uncertain location syntax:
        >>> variant = parse_hgvs("NC_000023.10:g.(?_32238146)_(32984039_?)del")
        >>> variant.description.location.start.start.is_unknown
        True
        >>> variant.description.location.end.end.coordinate is None
        True
    """

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
        """Return ``True`` when this coordinate has a known value."""
        return self.coordinate is not None

    @property
    def is_unknown(self) -> bool:
        """Return ``True`` when this coordinate is written as ``?``."""
        return self.coordinate is None

    @property
    def is_intronic(self) -> bool:
        """Return ``True`` for intronic coordinates.

        Examples:
            >>> from tinyhgvs import parse_hgvs
            >>> parse_hgvs("NM_004006.2:c.357+1G>A").description.location.start.is_intronic
            True
            >>> parse_hgvs("NM_001385026.1:c.-106+2T>A").description.location.start.is_intronic
            True
        """
        return self.offset not in (None, 0)

    @property
    def is_cds_start_anchored(self) -> bool:
        """Return ``True`` when the coordinate is relative to the CDS start."""
        return self.anchor is NucleotideAnchor.RELATIVE_CDS_START

    @property
    def is_cds_end_anchored(self) -> bool:
        """Return ``True`` when the coordinate is relative to the CDS end."""
        return self.anchor is NucleotideAnchor.RELATIVE_CDS_END

    @property
    def is_five_prime_utr(self) -> bool:
        """Return ``True`` for exonic positions in the 5' UTR.

        Examples:
            >>> from tinyhgvs import parse_hgvs
            >>> position = parse_hgvs("NM_007373.4:c.-123C>T").description.location.start
            >>> position.is_five_prime_utr
            True
            >>> position = parse_hgvs("NM_004006.2:c.76A>G").description.location.start
            >>> position.is_five_prime_utr
            False
        """
        return self.is_cds_start_anchored and self.offset == 0

    @property
    def is_three_prime_utr(self) -> bool:
        """Return ``True`` for exonic positions in the 3' UTR.

        Examples:
            >>> from tinyhgvs import parse_hgvs
            >>> position = parse_hgvs("NM_001272071.2:c.*1C>T").description.location.start
            >>> position.is_three_prime_utr
            True
            >>> position = parse_hgvs("NM_004006.2:c.76A>G").description.location.start
            >>> position.is_three_prime_utr
            False
        """
        return self.is_cds_end_anchored and self.offset == 0


@dataclass(frozen=True, slots=True)
class NucleotideNoChange:
    """Nucleotide no-change edit.

    Examples:
        >>> from tinyhgvs import NucleotideNoChange, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.2376=")
        >>> isinstance(variant.description, NucleotideNoChange)
        True
    """

    location: Location[NucleotideCoordinate]


@dataclass(frozen=True, slots=True)
class NucleotideSubstitution:
    """Nucleotide substitution.

    Examples:
        A reference base ``C`` is substituted by ``A``:
        >>> from tinyhgvs import NucleotideSubstitution, parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.33038255C>A")
        >>> isinstance(variant.description, NucleotideSubstitution)
        True
        >>> variant.description.reference
        'C'
        >>> variant.description.alternate
        'A'
    """

    location: Location[NucleotideCoordinate]
    reference: str
    alternate: str


@dataclass(frozen=True, slots=True)
class NucleotideDeletion:
    """Nucleotide deletion.

    Examples:
        >>> from tinyhgvs import NucleotideDeletion, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.5697del")
        >>> isinstance(variant.description, NucleotideDeletion)
        True
        >>> variant.description.location.start.coordinate
        5697
    """

    location: Location[NucleotideCoordinate]


@dataclass(frozen=True, slots=True)
class NucleotideDuplication:
    """Nucleotide duplication.

    Examples:
        >>> from tinyhgvs import NucleotideDuplication, parse_hgvs
        >>> variant = parse_hgvs("NC_000001.11:g.1234_2345dup")
        >>> isinstance(variant.description, NucleotideDuplication)
        True
        >>> variant.description.location.end.coordinate
        2345
    """

    location: Location[NucleotideCoordinate]


@dataclass(frozen=True, slots=True)
class NucleotideInversion:
    """Nucleotide inversion.

    Examples:
        >>> from tinyhgvs import NucleotideInversion, parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.32361330_32361333inv")
        >>> isinstance(variant.description, NucleotideInversion)
        True
    """

    location: Location[NucleotideCoordinate]


@dataclass(frozen=True, slots=True)
class NucleotideRepeat:
    """Top-level nucleotide repeat edit.

    Examples:
        A DNA repeat variant with a literal repeat unit:
        >>> from tinyhgvs import NucleotideRepeat, parse_hgvs
        >>> variant = parse_hgvs("NC_000014.8:g.123CAG[23]")
        >>> isinstance(variant.description, NucleotideRepeat)
        True
        >>> repeat = variant.description.sequence[0]
        >>> repeat.unit.value
        'CAG'
        >>> repeat.quantity.count
        23

        A RNA repeat variant can be composed of consecutive repeat blocks:
        >>> variant = parse_hgvs("NM_004006.3:r.456_499us[4]cag[9]gccag[3]")
        >>> len(variant.description.edit.sequence)
        3
        >>> variant.description.edit.sequence[2].quantity.count
        3
    """

    location: Location[NucleotideCoordinate]
    sequence: tuple[Repeat, ...]


@dataclass(frozen=True, slots=True)
class NucleotideInsertion:
    """Nucleotide insertion.

    Examples:
        Literal nucleotide insertion:
        >>> from tinyhgvs import LiteralSequence, NucleotideInsertion, parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.32862923_32862924insCCT")
        >>> isinstance(variant.description, NucleotideInsertion)
        True
        >>> item = variant.description.sequence[0]
        >>> isinstance(item, LiteralSequence)
        True
        >>> item.value
        'CCT'

        A composite insertion can mix literal and copied sequence:
        >>> variant = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")
        >>> len(variant.description.sequence)
        3
        >>> variant.description.sequence[0].value
        'T'
        >>> variant.description.sequence[1].is_from_same_reference
        True
        >>> variant.description.sequence[2].value
        'AGGG'

        Insertion of repeated unspecified bases uses the shared repeat model:
        >>> from tinyhgvs import UnknownRepeatUnit
        >>> variant = parse_hgvs("NC_000023.10:g.32717298_32717299insN[100]")
        >>> repeat = variant.description.sequence[0]
        >>> isinstance(repeat.unit, UnknownRepeatUnit)
        True
        >>> repeat.quantity.count
        100
    """

    location: Location[NucleotideCoordinate]
    sequence: tuple[NucleotideSequenceItem, ...]


@dataclass(frozen=True, slots=True)
class NucleotideDeletionInsertion:
    """Nucleotide deletion-insertion.

    Examples:
        A deleted interval is replaced by one literal sequence component:
        >>> from tinyhgvs import NucleotideDeletionInsertion, parse_hgvs
        >>> variant = parse_hgvs("LRG_199t1:c.850_901delinsTTCCTCGATGCCTG")
        >>> isinstance(variant.description, NucleotideDeletionInsertion)
        True
        >>> variant.description.sequence[0].value
        'TTCCTCGATGCCTG'

        Replacement sequence can be copied from the same reference:
        >>> variant = parse_hgvs("NC_000022.10:g.42522624_42522669delins42536337_42536382")
        >>> variant.description.sequence[0].location.start.coordinate
        42536337

        Replacement sequence can be a repeat item:
        >>> from tinyhgvs import UnknownRepeatUnit
        >>> variant = parse_hgvs("NM_004006.2:c.812_829delinsN[12]")
        >>> repeat = variant.description.sequence[0]
        >>> isinstance(repeat.unit, UnknownRepeatUnit)
        True
        >>> repeat.quantity.count
        12
    """

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
"""Tagged union for supported nucleotide edit models:

- [`NucleotideNoChange`][tinyhgvs.models.nucleotide.NucleotideNoChange]
- [`NucleotideSubstitution`][tinyhgvs.models.nucleotide.NucleotideSubstitution]
- [`NucleotideDeletion`][tinyhgvs.models.nucleotide.NucleotideDeletion]
- [`NucleotideDuplication`][tinyhgvs.models.nucleotide.NucleotideDuplication]
- [`NucleotideRepeat`][tinyhgvs.models.nucleotide.NucleotideRepeat]
- [`NucleotideInsertion`][tinyhgvs.models.nucleotide.NucleotideInsertion]
- [`NucleotideInversion`][tinyhgvs.models.nucleotide.NucleotideInversion]
- [`NucleotideDeletionInsertion`][tinyhgvs.models.nucleotide.NucleotideDeletionInsertion]
"""


@dataclass(frozen=True, slots=True)
class LiteralSequence:
    """Literal-base sequence component.

    Examples:
        >>> from tinyhgvs import LiteralSequence, parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.32862923_32862924insCCT")
        >>> item = variant.description.sequence[0]
        >>> isinstance(item, LiteralSequence)
        True
        >>> item.value
        'CCT'
    """

    value: str


@dataclass(frozen=True, slots=True)
class CopiedSequence:
    """Copied nucleotide sequence used in an insertion or deletion-insertion.

    Attributes:
        reference: Source reference when the copied sequence comes from a
            different accession. ``None`` means the same outer reference.
        coordinate_system: Source coordinate system when it differs from the
            outer variant. ``None`` means the same outer coordinate system.
        location: Location on the source reference.
        is_inverted: Whether the copied sequence is inserted in reverse
            orientation.

    Examples:
        A stretch of sequence from the same transcript is inserted in reverse
        orientation:
        >>> from tinyhgvs import CopiedSequence, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.849_850ins850_900inv")
        >>> item = variant.description.sequence[0]
        >>> isinstance(item, CopiedSequence)
        True
        >>> item.is_from_same_reference
        True
        >>> item.location.start.coordinate
        850
        >>> item.location.end.coordinate
        900
        >>> item.is_inverted
        True

        A copied sequence can also come from another chromosome:
        >>> variant = parse_hgvs(
        ...     "NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]"
        ... )
        >>> item = variant.description.sequence[0]
        >>> item.reference.primary.id
        'NC_000022.10'
        >>> item.coordinate_system
        <CoordinateSystem.GENOMIC: 'g'>
    """

    reference: ReferenceSpec | None
    coordinate_system: CoordinateSystem | None
    location: KnownLocation[NucleotideCoordinate]
    is_inverted: bool

    @property
    def is_from_same_reference(self) -> bool:
        return self.reference is None and self.coordinate_system is None


NucleotideSequenceItem: TypeAlias = LiteralSequence | Repeat | CopiedSequence
"""Tagged union for supported inserted or replacement nucleotide components:

- [`LiteralSequence`][tinyhgvs.models.nucleotide.LiteralSequence]
- [`Repeat`][tinyhgvs.models.repeat.Repeat]
- [`CopiedSequence`][tinyhgvs.models.nucleotide.CopiedSequence]
"""


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
