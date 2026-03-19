"""Shared Python model types for parsed HGVS variants.

This module groups the model pieces used by both nucleotide and protein
variants: reference identifiers, coordinate-system labels, and generic
intervals over a reference.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Generic, TypeVar


class CoordinateSystem(str, Enum):
    """Supported HGVS coordinate types.

    The coordinate system tells users what kind of biological reference frame
    is being used by the parsed variant.

    Attributes:
        GENOMIC: Genomic DNA coordinates written as ``g.``.
        CIRCULAR_GENOMIC: Circular genomic DNA coordinates written as ``o.``.
        MITOCHONDRIAL: Mitochondrial DNA coordinates written as ``m.``.
        CODING_DNA: Coding DNA coordinates written as ``c.``.
        NON_CODING_DNA: Non-coding DNA coordinates written as ``n.``.
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
    CIRCULAR_GENOMIC = "o"
    MITOCHONDRIAL = "m"
    CODING_DNA = "c"
    NON_CODING_DNA = "n"
    RNA = "r"
    PROTEIN = "p"


@dataclass(frozen=True, slots=True)
class Accession:
    """Sequence accession with optional version.

    Attributes:
        id: Accession string exactly as it appears in the HGVS expression.
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


PositionT = TypeVar("PositionT")


@dataclass(frozen=True, slots=True)
class Interval(Generic[PositionT]):
    """Inclusive interval used for nucleotide and protein locations.

    When ``end`` is omitted, the interval represents a single coordinate.

    Attributes:
        start: Start position in the range.
        end: Optional inclusive end position.

    Examples:
        A single-position nucleotide substitution has no end coordinate:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> variant_location = variant.description.location
        >>> variant_location.start.coordinate
        357
        >>> variant_location.end is None
        True

        A protein deletion spanning multiple residues has both start and end:
        >>> variant = parse_hgvs("NP_003997.2:p.Lys23_Val25del")
        >>> variant_location = variant.description.effect.location
        >>> variant_location.start.residue
        'Lys'
        >>> variant_location.end.residue
        'Val'
    """

    start: PositionT
    end: PositionT | None = None


__all__ = [
    "Accession",
    "CoordinateSystem",
    "Interval",
    "PositionT",
    "ReferenceSpec",
]
