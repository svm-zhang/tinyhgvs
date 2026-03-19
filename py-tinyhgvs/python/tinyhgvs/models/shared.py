"""Shared Python model types for parsed HGVS variants.

This module groups the model pieces used by both nucleotide and protein
variants: reference identifiers, coordinate-system labels, and generic ranges
over a reference.

Type Aliases:
    SequenceSource: Tagged union for all supported sequence edit source models.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Generic, Literal, TypeAlias, TypeVar


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


class SequenceKind(str, Enum):
    """Accession families inferred from sequence identifiers.

    Attributes:
        REFSEQ_CHROMOSOME: RefSeq chromosome accession such as ``NC_000023.11``.
        REFSEQ_CONTIG: RefSeq contig accession such as ``NT_`` or ``NW_``.
        REFSEQ_GENE_REGION: RefSeq gene-region accession such as ``NG_``.
        REFSEQ_CODING_TRANSCRIPT: RefSeq coding transcript accession such as ``NM_``.
        REFSEQ_NONCODING_TRANSCRIPT: RefSeq non-coding transcript accession such as ``NR_``.
        REFSEQ_PROTEIN: RefSeq protein accession such as ``NP_``.
        ENSEMBL_GENE: Ensembl gene accession such as ``ENSG...``.
        ENSEMBL_TRANSCRIPT: Ensembl transcript accession such as ``ENST...``.
        ENSEMBL_PROTEIN: Ensembl protein accession such as ``ENSP...``.
        LRG_GENE_REGION: LRG genomic accession.
        LRG_TRANSCRIPT: LRG transcript accession.
        LRG_PROTEIN: LRG protein accession.
        UNKNOWN: Accession family not currently recognized by the parser.

    Examples:
        A RefSeq coding transcript accession:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.5657_5660inv")
        >>> variant.reference.primary.kind
        <SequenceKind.REFSEQ_CODING_TRANSCRIPT: 'refseq_coding_transcript'>

        A LRG transcript accession:
        >>> variant = parse_hgvs("LRG_199t1:c.850_901delinsTTCCTCGATGCCTG")
        >>> variant.reference.primary.kind
        <SequenceKind.LRG_TRANSCRIPT: 'lrg_transcript'>

        An Ensembl gene accession:
        >>> variant = parse_hgvs("ENSG00000160190.9(ENST00000352133.2):c.1521+898G>A")
        >>> variant.reference.primary.kind
        <SequenceKind.ENSEMBL_GENE: 'ensembl_gene'>
    """

    REFSEQ_CHROMOSOME = "refseq_chromosome"
    REFSEQ_CONTIG = "refseq_contig"
    REFSEQ_GENE_REGION = "refseq_gene_region"
    REFSEQ_CODING_TRANSCRIPT = "refseq_coding_transcript"
    REFSEQ_NONCODING_TRANSCRIPT = "refseq_noncoding_transcript"
    REFSEQ_PROTEIN = "refseq_protein"
    ENSEMBL_GENE = "ensembl_gene"
    ENSEMBL_TRANSCRIPT = "ensembl_transcript"
    ENSEMBL_PROTEIN = "ensembl_protein"
    LRG_GENE_REGION = "lrg_gene_region"
    LRG_TRANSCRIPT = "lrg_transcript"
    LRG_PROTEIN = "lrg_protein"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class SequenceId:
    """Sequence identifier with optional version and inferred accession family.

    Attributes:
        raw: Raw accession string.
        version: Numeric suffix when a versioned accession is present.
        kind: Inferred accession family.

    Examples:
        A protein accession with version:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NP_003997.2:p.Val7del")
        >>> variant.reference.primary.raw
        'NP_003997.2'
        >>> variant.reference.primary.version
        2
        >>> variant.reference.primary.kind
        <SequenceKind.REFSEQ_PROTEIN: 'refseq_protein'>

        A gene accession with transcript context:
        >>> variant = parse_hgvs("ENSG00000160190.9(ENST00000352133.2):c.1521+898G>A")
        >>> variant.reference.primary.raw
        'ENSG00000160190.9'
        >>> variant.reference.primary.version
        9
        >>> variant.reference.primary.kind
        <SequenceKind.ENSEMBL_GENE: 'ensembl_gene'>
    """

    raw: str
    version: int | None
    kind: SequenceKind


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
        >>> variant.reference.primary.raw
        'NC_000023.10'
        >>> variant.reference.context is None
        True

        A coding variant described on a genomic reference with transcript
        context:
        >>> variant = parse_hgvs("NC_000023.11(NM_004006.2):c.3921dup")
        >>> variant.reference.primary.raw
        'NC_000023.11'
        >>> variant.reference.context.raw
        'NM_004006.2'
    """

    primary: SequenceId
    context: SequenceId | None


PositionT = TypeVar("PositionT")


@dataclass(frozen=True, slots=True)
class Range(Generic[PositionT]):
    """Inclusive range model for positions that locate a variant.

    ``Range`` is used for both nucleotide and protein locations. When ``end``
    is omitted, the range represents a single position.

    Attributes:
        start: Start position in the range.
        end: Optional inclusive end position.

    Examples:
        A single-position nucleotide substitution has no end coordinate:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> variant_location = variant.description.location
        >>> variant_location.start.position
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


@dataclass(frozen=True, slots=True)
class CurrentReferenceSource:
    """Model describing a sequence segment copied from the same reference.

    The inserted or replacement sequence is copied from the same primary
    reference that the HGVS variant itself uses.

    Attributes:
        kind: Marks the source as the current reference.

    Examples:
        A segment of sequence from the same reference ``NM_004006.2`` between
        positions 858 and 895 is inserted. The edit has one single component.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.849_850ins858_895")
        >>> variant_edit_component = variant.description.edit.sequence.components[0]
        >>> variant_edit_component.segment.source
        CurrentReferenceSource(kind='current_reference')
    """

    kind: Literal["current_reference"] = field(
        init=False,
        default="current_reference",
    )


@dataclass(frozen=True, slots=True)
class OtherReferenceSource:
    """Model describing a sequence segment copied from another reference.

    Attributes:
        reference: Reference sequence from which the segment is copied.
        coordinate_system: Coordinate system used on that external reference.
        kind: Marks the source as a different reference.

    Examples:
        A replacement segment is copied from chromosome 22 and inserted into a
        variant described on chromosome 12. The edit has one segment component.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs(
        ...     "NC_000012.11:g.6128892_6128954delins[NC_000022.10:g.17179029_17179091]"
        ... )
        >>> variant_edit_component = variant.description.edit.sequence.components[0]
        >>> variant_edit_source = variant_edit_component.segment.source
        >>> variant_edit_source.reference.primary.raw
        'NC_000022.10'
        >>> variant_edit_source.reference.primary.raw != variant.reference.primary.raw
        True
    """

    reference: ReferenceSpec
    coordinate_system: CoordinateSystem
    kind: Literal["other_reference"] = field(
        init=False,
        default="other_reference",
    )


SequenceSource: TypeAlias = CurrentReferenceSource | OtherReferenceSource
"""Tagged union for all supported sequence edit source models:

- [`CurrentReferenceSource`][tinyhgvs.models.shared.CurrentReferenceSource]
- [`OtherReferenceSource`][tinyhgvs.models.shared.OtherReferenceSource]
"""


__all__ = [
    "CoordinateSystem",
    "CurrentReferenceSource",
    "OtherReferenceSource",
    "PositionT",
    "Range",
    "ReferenceSpec",
    "SequenceId",
    "SequenceKind",
    "SequenceSource",
]
