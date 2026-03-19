"""Nucleotide-focused Python model types for parsed HGVS variants.

Type Aliases:
    NucleotideSequenceComponent: Tagged union for supported inserted or
        replacement nucleotide sequence component models.
    NucleotideEdit: Tagged union for supported nucleotide edit models.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Literal, TypeAlias

from .shared import (
    CurrentReferenceSource,
    Range,
    SequenceSource,
)


class NucleotidePositionAnchor(str, Enum):
    """HGVS reference point used to interpret a nucleotide position.

    Attributes:
        COORDINATE: Position is interpreted directly on the reference sequence.
        CDS_START: Position is interpreted relative to the coding DNA start
            site.
        CDS_END: Position is interpreted relative to the coding DNA end site.

    Examples:
        An intronic splice-site substitution uses direct coordinates:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> variant.description.location.start.anchor
        <NucleotidePositionAnchor.COORDINATE: 'coordinate'>

        A 5' UTR substitution is anchored to the CDS start:
        >>> variant = parse_hgvs("NM_007373.4:c.-1C>T")
        >>> variant.description.location.start.anchor
        <NucleotidePositionAnchor.CDS_START: 'cds_start'>

        A 3' UTR substitution is anchored to the CDS end:
        >>> variant = parse_hgvs("NM_001272071.2:c.*1C>T")
        >>> variant.description.location.start.anchor
        <NucleotidePositionAnchor.CDS_END: 'cds_end'>
    """

    COORDINATE = "coordinate"
    CDS_START = "cds_start"
    CDS_END = "cds_end"


@dataclass(frozen=True, slots=True)
class NucleotidePosition:
    """Nucleotide position with anchor and signed offset semantics.

    Attributes:
        anchor: Reference point used to interpret the position.
        position: Concrete coordinate when it is known. For 5' UTR positions
            such as ``c.-81``, this is normalized to ``0``. For 3' UTR
            positions such as ``c.*24``, the terminal CDS coordinate cannot be
            inferred from the string alone, so this is ``None``.
        offset: Signed offset from ``position`` or the anchor. Positive values
            move downstream, and negative values move upstream.

    Examples:
        Duplication crossing an exon/intron border:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.11(NM_004006.2):c.260_264+48dup")
        >>> variant.description.location.start.position
        260
        >>> variant.description.location.end.position
        264
        >>> variant.description.location.end.offset
        48

        Upstream intronic substitution:
        >>> variant = parse_hgvs("NG_012232.1(NM_004006.2):c.264-2A>G")
        >>> variant.description.location.start.position
        264
        >>> variant.description.location.start.offset
        -2

        5' UTR substitution:
        >>> variant = parse_hgvs("NM_007373.4:c.-1C>T")
        >>> variant.description.location.start.anchor
        <NucleotidePositionAnchor.CDS_START: 'cds_start'>
        >>> variant.description.location.start.position
        0
        >>> variant.description.location.start.offset
        -1

        3' UTR substitution:
        >>> variant = parse_hgvs("NM_001272071.2:c.*1C>T")
        >>> variant.description.location.start.anchor
        <NucleotidePositionAnchor.CDS_END: 'cds_end'>
        >>> variant.description.location.start.position is None
        True
        >>> variant.description.location.start.offset
        1
    """

    anchor: NucleotidePositionAnchor
    position: int | None
    offset: int = 0

    @property
    def is_intronic(self) -> bool:
        """Return ``True`` for coordinate-anchored positions with an offset.

        Examples:
            >>> from tinyhgvs import parse_hgvs
            >>> parse_hgvs("NM_004006.2:c.357+1G>A").description.location.start.is_intronic
            True
            >>> parse_hgvs("NG_012232.1(NM_004006.2):c.264-2A>G").description.location.start.is_intronic
            True
        """
        return (
            self.anchor is NucleotidePositionAnchor.COORDINATE
            and self.offset != 0
        )

    @property
    def _is_cds_start_relative(self) -> bool:
        return self.anchor is NucleotidePositionAnchor.CDS_START

    @property
    def _is_cds_end_relative(self) -> bool:
        return self.anchor is NucleotidePositionAnchor.CDS_END

    @property
    def is_five_prime_utr(self) -> bool:
        """Return ``True`` for positions in the 5' UTR.

        Examples:
            >>> from tinyhgvs import parse_hgvs
            >>> position = parse_hgvs("NM_007373.4:c.-123C>T").description.location.start
            >>> position.is_five_prime_utr
            True
            >>> position.is_three_prime_utr
            False
        """
        return self._is_cds_start_relative

    @property
    def is_three_prime_utr(self) -> bool:
        """Return ``True`` for positions in the 3' UTR.

        Examples:
            >>> from tinyhgvs import parse_hgvs
            >>> position = parse_hgvs("NM_001272071.2:c.*1C>T").description.location.start
            >>> position.is_three_prime_utr
            True
            >>> position.is_five_prime_utr
            False
        """
        return self._is_cds_end_relative


@dataclass(frozen=True, slots=True)
class NucleotideSequenceSegment:
    """Segment-type inserted or replacement sequence component.

    Attributes:
        source: Reference source from which the segment is copied.
        location: Inclusive location on that source reference.
        is_inverted: Whether the copied segment is inverted before insertion.

    Examples:
        A segment of sequence from the same reference is inserted after being
        reversed in orientation. The edit has one segment component.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.849_850ins850_900inv")
        >>> variant_edit_component = variant.description.edit.sequence.components[0]
        >>> variant_edit_segment = variant_edit_component.segment
        >>> variant_edit_segment.source
        CurrentReferenceSource(kind='current_reference')
        >>> variant_edit_segment.location.start.position
        850
        >>> variant_edit_segment.location.end.position
        900
        >>> variant_edit_segment.is_inverted
        True
    """

    source: SequenceSource
    location: Range[NucleotidePosition]
    is_inverted: bool


@dataclass(frozen=True, slots=True)
class LiteralSequenceComponent:
    """Model describing literal-base-type sequence edit component.

    Attributes:
        value: Nucleotide bases.
        kind: Marks the sequence edit component as literal nucleotide bases.

    Examples:
        A literal insertion of three nucleotides. The edit has one single
        literal component.
        >>> from tinyhgvs import parse_hgvs, LiteralSequenceComponent
        >>> variant = parse_hgvs("NC_000023.10:g.32862923_32862924insCCT")
        >>> variant_edit_component = variant.description.edit.sequence.components[0]
        >>> isinstance(variant_edit_component, LiteralSequenceComponent)
        True
        >>> variant_edit_component.value
        'CCT'
        >>> variant_edit_component.kind
        'literal'
    """

    value: str
    kind: Literal["literal"] = field(init=False, default="literal")


@dataclass(frozen=True, slots=True)
class RepeatSequenceComponent:
    """Model describing repeat-type sequence edit component.

    Attributes:
        unit: Repeat unit of nucleotide bases.
        count: Number of units being repeated.
        kind: Marks the sequence edit component as repeat-type.

    Examples:
        The insertion contains one repeated component of unspecified nucleotides.
        >>> from tinyhgvs import parse_hgvs, RepeatSequenceComponent
        >>> variant = parse_hgvs("NC_000023.10:g.32717298_32717299insN[100]")
        >>> variant_edit_component = variant.description.edit.sequence.components[0]
        >>> isinstance(variant_edit_component, RepeatSequenceComponent)
        True
        >>> variant_edit_component.unit
        'N'
        >>> variant_edit_component.count
        100
        >>> variant_edit_component.kind
        'repeat'
    """

    unit: str
    count: int
    kind: Literal["repeat"] = field(init=False, default="repeat")


@dataclass(frozen=True, slots=True)
class SegmentSequenceComponent:
    """Model describing segment-type sequence edit component.

    Attributes:
        segment: Segment modeled by :class:`NucleotideSequenceSegment`.
        kind: Marks the sequence edit component as segment-type.

    Examples:
        The deleted reference sequence is replaced by one copied segment from
        the same reference. The edit has one segment component.
        >>> from tinyhgvs import parse_hgvs, NucleotideSequenceSegment
        >>> variant = parse_hgvs("NC_000022.10:g.42522624_42522669delins42536337_42536382")
        >>> variant_edit_component = variant.description.edit.sequence.components[0]
        >>> isinstance(variant_edit_component.segment, NucleotideSequenceSegment)
        True
        >>> variant_edit_component.kind
        'segment'
    """

    segment: NucleotideSequenceSegment
    kind: Literal["segment"] = field(init=False, default="segment")


NucleotideSequenceComponent: TypeAlias = (
    LiteralSequenceComponent
    | RepeatSequenceComponent
    | SegmentSequenceComponent
)
"""Tagged union for supported inserted or replacement nucleotide components: 

- [`LiteralSequenceComponent`][tinyhgvs.models.nucleotide.LiteralSequenceComponent]
- [`RepeatSequenceComponent`][tinyhgvs.models.nucleotide.RepeatSequenceComponent]
- [`SegmentSequenceComponent`][tinyhgvs.models.nucleotide.SegmentSequenceComponent]
"""


@dataclass(frozen=True, slots=True)
class NucleotideSequence:
    """Model describing inserted and replacement nucleotide sequence.

    Attributes:
        components: Tuple of sequence components modeled by
            :class:`NucleotideSequenceComponent`.

    Examples:
        A deletion-insertion with one literal replacement component:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("LRG_199t1:c.850_901delinsTTCCTCGATGCCTG")
        >>> len(variant.description.edit.sequence.components)
        1

        An insertion with a literal component, a copied segment, and another
        literal component:
        >>> variant = parse_hgvs("LRG_199t1:c.419_420ins[T;450_470;AGGG]")
        >>> len(variant.description.edit.sequence.components)
        3
    """

    components: tuple[NucleotideSequenceComponent, ...]


@dataclass(frozen=True, slots=True)
class NucleotideNoChangeEdit:
    """Model describing no nucleotide change.

    Attributes:
        kind: Edit kind.

    Examples:
        A coding DNA position is explicitly reported as unchanged.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.123=")
        >>> variant.description.edit.kind
        'no_change'
    """

    kind: Literal["no_change"] = field(init=False, default="no_change")


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
class NucleotideDeletionEdit:
    """Model describing nucleotide deletion.

    Attributes:
        kind: Edit kind.

    Examples:
        A nucleotide interval spanning intronic sequence is deleted.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.4072-1234_5155-246del")
        >>> variant.description.edit.kind
        'deletion'
    """

    kind: Literal["deletion"] = field(init=False, default="deletion")


@dataclass(frozen=True, slots=True)
class NucleotideDuplicationEdit:
    """Model describing nucleotide duplication.

    Attributes:
        kind: Edit kind.

    Examples:
        A genomic interval is duplicated.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000001.11:g.1234_2345dup")
        >>> variant.description.edit.kind
        'duplication'
    """

    kind: Literal["duplication"] = field(init=False, default="duplication")


@dataclass(frozen=True, slots=True)
class NucleotideInsertionEdit:
    """Model describing nucleotide insertion.

    Attributes:
        sequence: Inserted nucleotide sequence modeled by
            :class:`NucleotideSequence`.
        kind: Edit kind.

    Examples:
        Literal nucleotide insertion:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.32862923_32862924insCCT")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.sequence.components[0].value
        'CCT'

        A composite insertion of a literal nucleotide followed by a segment
        copied from the same transcript reference. The edit has two components.
        >>> variant = parse_hgvs("NM_004006.2:c.419_420ins[T;401_419]")
        >>> variant_edit = variant.description.edit
        >>> len(variant_edit.sequence.components)
        2
        >>> variant_edit.sequence.components[0].kind
        'literal'
        >>> variant_edit.sequence.components[1].kind
        'segment'

        Insertion of unspecified repeated bases:
        >>> variant = parse_hgvs("NC_000023.10:g.32717298_32717299insN[100]")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.sequence.components[0].unit
        'N'
        >>> variant_edit.sequence.components[0].count
        100
    """

    sequence: NucleotideSequence
    kind: Literal["insertion"] = field(init=False, default="insertion")


@dataclass(frozen=True, slots=True)
class NucleotideInversionEdit:
    """Model describing nucleotide inversion.

    Attributes:
        kind: Edit kind.

    Examples:
        A short nucleotide interval is inverted in place.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.5657_5660inv")
        >>> variant.description.edit.kind
        'inversion'
    """

    kind: Literal["inversion"] = field(init=False, default="inversion")


@dataclass(frozen=True, slots=True)
class NucleotideDeletionInsertionEdit:
    """Model describing nucleotide deletion-insertion.

    Attributes:
        sequence: Replacement nucleotide sequence modeled by
            :class:`NucleotideSequence`.
        kind: Edit kind.

    Examples:
        A deleted interval is replaced by one literal sequence component.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("LRG_199t1:c.850_901delinsTTCCTCGATGCCTG")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.sequence.components[0].value
        'TTCCTCGATGCCTG'

        A deleted interval is replaced by a copied segment from the same
        reference. The replacement has one segment component.
        >>> variant = parse_hgvs("NC_000022.10:g.42522624_42522669delins42536337_42536382")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.sequence.components[0].segment.location.start.position
        42536337

        A deleted interval is replaced by repeated unspecified bases.
        >>> variant = parse_hgvs("NM_004006.2:c.812_829delinsN[12]")
        >>> variant_edit = variant.description.edit
        >>> variant_edit.sequence.components[0].unit
        'N'
        >>> variant_edit.sequence.components[0].count
        12
    """

    sequence: NucleotideSequence
    kind: Literal["deletion_insertion"] = field(
        init=False,
        default="deletion_insertion",
    )


NucleotideEdit: TypeAlias = (
    NucleotideNoChangeEdit
    | NucleotideSubstitutionEdit
    | NucleotideDeletionEdit
    | NucleotideDuplicationEdit
    | NucleotideInsertionEdit
    | NucleotideInversionEdit
    | NucleotideDeletionInsertionEdit
)
"""Tagged union for supported nucleotide edit models:

- [`NucleotideNoChangeEdit`][tinyhgvs.models.nucleotide.NucleotideNoChangeEdit]
- [`NucleotideSubstitutionEdit`][tinyhgvs.models.nucleotide.NucleotideSubstitutionEdit]
- [`NucleotideDeletionEdit`][tinyhgvs.models.nucleotide.NucleotideDeletionEdit]
- [`NucleotideDuplicationEdit`][tinyhgvs.models.nucleotide.NucleotideDuplicationEdit]
- [`NucleotideInsertionEdit`][tinyhgvs.models.nucleotide.NucleotideInsertionEdit]
- [`NucleotideInversionEdit`][tinyhgvs.models.nucleotide.NucleotideInversionEdit]
- [`NucleotideDeletionInsertionEdit`][tinyhgvs.models.nucleotide.NucleotideDeletionInsertionEdit]
"""


@dataclass(frozen=True, slots=True)
class NucleotideVariant:
    """Model describing a nucleotide-level variant.

    Attributes:
        location: Inclusive nucleotide interval where the edit is applied.
        edit: Nucleotide edit applied at that interval.

    Examples:
        A splice-site substitution is represented by a nucleotide location and
        a nucleotide substitution edit.
        >>> from tinyhgvs import NucleotideSubstitutionEdit, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> isinstance(variant.description.edit, NucleotideSubstitutionEdit)
        True
        >>> variant_description = variant.description
        >>> variant_description.location.start.position
        357
        >>> variant_description.location.start.offset
        1
    """

    location: Range[NucleotidePosition]
    edit: NucleotideEdit


__all__ = [
    "CoordinateSystem",
    "CurrentReferenceSource",
    "NucleotideDeletionEdit",
    "NucleotideDeletionInsertionEdit",
    "NucleotideDuplicationEdit",
    "NucleotideEdit",
    "NucleotideInsertionEdit",
    "NucleotideInversionEdit",
    "NucleotideNoChangeEdit",
    "NucleotidePosition",
    "NucleotidePositionAnchor",
    "NucleotideSequence",
    "NucleotideSequenceComponent",
    "NucleotideSequenceSegment",
    "NucleotideSubstitutionEdit",
    "NucleotideVariant",
    "RepeatSequenceComponent",
    "SegmentSequenceComponent",
]
