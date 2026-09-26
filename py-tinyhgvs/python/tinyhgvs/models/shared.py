"""Shared Python model types for parsed HGVS variants.

This module groups the model pieces used by both nucleotide and protein
variants: reference identifiers, coordinate-system labels, generic intervals,
and allele container types.
"""

from __future__ import annotations, with_statement

from dataclasses import dataclass
from enum import Enum
from typing import Generic, TypeAlias, TypeVar


class CoordinateSystem(str, Enum):
    """Supported HGVS coordinate types.

    The coordinate system tells users what kind of biological reference frame
    is being used by the parsed variant.

    Attributes:
        GENOMIC: Genomic DNA coordinates written as ``g.``.
        CODING_DNA: Coding DNA coordinates written as ``c.``.
        RNA: RNA coordinates written as ``r.``.
        PROTEIN: Protein coordinates written as ``p.``.

    Examples:
        Genomic DNA variant:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.11:g.33038255C>A")
        >>> variant.coordinate_system
        <CoordinateSystem.GENOMIC: 'g'>

        Coding DNA variant:
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> variant.coordinate_system
        <CoordinateSystem.CODING_DNA: 'c'>

        Protein variant:
        >>> variant = parse_hgvs("NP_003997.1:p.Trp24Ter")
        >>> variant.coordinate_system
        <CoordinateSystem.PROTEIN: 'p'>
    """

    GENOMIC = "g"
    CODING_DNA = "c"
    RNA = "r"
    PROTEIN = "p"


_PositionT = TypeVar("_PositionT")


@dataclass(frozen=True, slots=True)
class Accession:
    """Sequence accession with optional version.

    Attributes:
        id: Accession string as it appears in the HGVS expression.
        version: Parsed version suffix when one is present.

    Examples:
        A RefSeq protein accession with version:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NP_003997.2:p.Val7del")
        >>> variant.reference.primary.id
        'NP_003997.2'
        >>> variant.reference.primary.version
        2

        A transcript accession without an explicit contextual reference:
        >>> variant = parse_hgvs("NM_007373.4:c.-1C>T")
        >>> variant.reference.primary.id
        'NM_007373.4'
        >>> variant.reference.primary.version
        4

        A genomic accession with transcript context:
        >>> variant = parse_hgvs("ENSG00000160190.9(ENST00000352133.2):c.1521+898G>A")
        >>> variant.reference.primary.id
        'ENSG00000160190.9'
        >>> variant.reference.primary.version
        9
    """

    id: str
    version: int | None


@dataclass(frozen=True, slots=True)
class ReferenceSpec:
    """Reference sequence field preceding the ``:`` in a HGVS string.

    Attributes:
        primary: Primary accession being described.
        context: Optional contextual accession, commonly a transcript nested
            inside a genomic reference.

    Examples:
        A variant described directly on one reference sequence:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.33038255C>A")
        >>> variant.reference.primary.id
        'NC_000023.10'
        >>> variant.reference.context is None
        True

        A coding variant described on a genomic reference with transcript
        context:
        >>> variant = parse_hgvs("NC_000023.11(NM_004006.2):c.3921dup")
        >>> variant.reference.primary.id
        'NC_000023.11'
        >>> variant.reference.context.id
        'NM_004006.2'
    """

    primary: Accession
    context: Accession | None


@dataclass(frozen=True, slots=True)
class KnownLocation(Generic[_PositionT]):
    """A known HGVS location.

    A missing end represents a single position.

    Examples:
        A single-position coding DNA substitution has no end coordinate:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> location = variant.description.location
        >>> location.start.coordinate
        357
        >>> location.end is None
        True

        A protein deletion spanning multiple residues has both start and end:
        >>> variant = parse_hgvs("NP_003997.2:p.Lys23_Val25del")
        >>> location = variant.description.edit.location
        >>> location.start.residue
        'Lys'
        >>> location.end.residue
        'Val'
    """

    start: _PositionT
    end: _PositionT | None = None

    @property
    def is_position(self) -> bool:
        return self.end is None

    @property
    def is_interval(self) -> bool:
        return self.end is not None


@dataclass(frozen=True, slots=True)
class PossibleRange(Generic[_PositionT]):
    """Range within which one uncertain location boundary may lie.

    Examples:
        A genomic deletion can have an uncertain left boundary:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.(?_32238146)_(32984039_?)del")
        >>> left = variant.description.location.start
        >>> left.start.is_unknown
        True
        >>> left.end.coordinate
        32238146
    """

    start: _PositionT
    end: _PositionT | None = None


@dataclass(frozen=True, slots=True)
class UncertainLocation(Generic[_PositionT]):
    """An HGVS location with uncertain positional boundaries.

    Examples:
        A substitution written in parentheses has an uncertain location:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.(33038277_33038278)C>T")
        >>> location = variant.description.location
        >>> location.is_interval
        False
        >>> location.start.start.coordinate
        33038277
        >>> location.start.end.coordinate
        33038278
        >>> location.end is None
        True

        A deletion can have uncertain left and right boundaries:
        >>> variant = parse_hgvs("NC_000023.10:g.(?_32238146)_(32984039_?)del")
        >>> location = variant.description.location
        >>> location.is_interval
        True
        >>> location.start.start.is_unknown
        True
        >>> location.start.end.coordinate
        32238146
        >>> location.end.start.coordinate
        32984039
        >>> location.end.end.is_unknown
        True

        Protein locations can also be uncertain:
        >>> variant = parse_hgvs("NP_003997.1:p.(Ala123_Pro131)Ter")
        >>> location = variant.description.edit.location
        >>> location.start.start.residue
        'Ala'
        >>> location.start.end.residue
        'Pro'
    """

    start: PossibleRange[_PositionT]
    end: PossibleRange[_PositionT] | None = None

    @property
    def is_position(self) -> bool:
        return self.end is None

    @property
    def is_interval(self) -> bool:
        return self.end is not None


Location: TypeAlias = KnownLocation[_PositionT] | UncertainLocation[_PositionT]
"""Tagged union for supported location models:

- [`KnownLocation`][tinyhgvs.models.shared.KnownLocation]
- [`UncertainLocation`][tinyhgvs.models.shared.UncertainLocation]
"""


class OutcomeCertainty(str, Enum):
    CERTAIN = "certain"
    PREDICTED = "predicted"


__all__ = [
    "Accession",
    "CoordinateSystem",
    "KnownLocation",
    "Location",
    "ReferenceSpec",
    "OutcomeCertainty",
    "UncertainLocation",
    "PossibleRange",
]
