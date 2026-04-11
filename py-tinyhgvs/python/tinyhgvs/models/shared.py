"""Shared Python model types for parsed HGVS variants.

This module groups the model pieces used by both nucleotide and protein
variants: reference identifiers, coordinate-system labels, generic intervals,
and allele container types.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Generic, Iterator, TypeVar


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
VariantT = TypeVar("VariantT")


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


class AllelePhase(str, Enum):
    """Model describing phase between/among alleles.

    Attributes:
        TRANS: Alleles are *In-trans* phase.
        UNCERTAIN: Phase between/among alleles is uncertain.

    Examples:
        *In-trans* alleles:
        >>> from tinyhgvs import parse_hgvs
        >>> variants = parse_hgvs("NM_004006.2:c.[2376G>C];[3103del]")
        >>> variants.description.phase
        <AllelePhase.TRANS: 'trans'>

        Uncertain phase:
        >>> variant = parse_hgvs("NC_000001.11:g.123G>A(;)345del")
        >>> variant.description.phase
        <AllelePhase.UNCERTAIN: 'uncertain'>

        *In-trans* protein alleles:
        >>> variants = parse_hgvs("NP_003997.1:p.[Ser68Arg];[Ser68=]")
        >>> variants.description.phase
        <AllelePhase.TRANS: 'trans'>

        Protein alleles with uncertain phase:
        >>> variant = parse_hgvs("NP_003997.1:p.(Ser73Arg)(;)(Asn103del)")
        >>> variant.description.phase
        <AllelePhase.UNCERTAIN: 'uncertain'>
    """

    TRANS = "trans"
    UNCERTAIN = "uncertain"


@dataclass(frozen=True, slots=True)
class Allele(Generic[VariantT]):
    """One allele carrying one or more variants *in cis*.

    Attributes:
        variants: Variants described on the same allele and therefore treated
            as occurring together *in cis*.

    Examples:
        A nucleotide allele carrying multiple variants:

        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000001.11:g.[123G>A;345del]")
        >>> len(variant.description.allele_one.variants)
        2

        A protein allele carrying multiple variants together:

        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.[Ser68Arg;Asn594del]")
        >>> len(variant.description.allele_one.variants)
        2
    """

    variants: tuple[VariantT, ...]

    def __iter__(self) -> Iterator[VariantT]:
        """Return an iterator over variants carried by an allele in order.

        Returns:
            (Iterator[VariantT]): Iterator over variants carried by this allele,
                in the order they appear in the written description.
        Examples:
            >>> from tinyhgvs import parse_hgvs
            >>> allele = parse_hgvs(
            ...     "NP_003997.1:p.[Ser68Arg;Asn594del]"
            ... ).description.allele_one
            >>> len(tuple(allele))
            2

            >>> from tinyhgvs import parse_hgvs
            >>> allele = parse_hgvs(
            ...     "NC_000001.11:g.[123G>A;345del]"
            ... ).description.allele_one
            >>> len(tuple(allele))
            2
        """
        return iter(self.variants)


@dataclass(frozen=True, slots=True)
class AlleleVariant(Generic[VariantT]):
    """Structured representation of HGVS allele variant syntax.

    An allele variant may describe:

    - a single allele carrying one or more variants.
    - two alleles with an explicit phase relationship.
    - additional alleles whose relation to the established allele state is
      uncertain.

    Attributes:
        allele_one: First allele in the allele-variant description.
        allele_two: Second allele, if present.
        phase: Phase relation between ``allele_one`` and ``allele_two``. This
            is ``None`` when only one allele is described.
        alleles_unphased: Additional alleles in uncertain relation to the
            allele state established by ``allele_one`` and ``allele_two``.

    Examples:
        Variants *in cis* on a single allele:

        >>> from tinyhgvs import parse_hgvs
        >>> desc = parse_hgvs("NC_000023.10:g.[30683643A>G;33038273T>G]").description
        >>> len(tuple(desc))
        1
        >>> desc.allele_two is None
        True
        >>> desc.phase is None
        True
        >>> len(desc.allele_one.variants)
        2

        Two nucleotide alleles *in trans*:

        >>> desc = parse_hgvs("NM_004006.2:c.[2376G>C];[3103del]").description
        >>> desc.allele_two is not None
        True
        >>> desc.phase
        <AllelePhase.TRANS: 'trans'>
        >>> len(desc.allele_one.variants)
        1
        >>> len(desc.allele_two.variants)
        1

        Additional alleles with uncertain phase:

        >>> desc = parse_hgvs(
        ...     "NM_004006.2:c.[296T>G;476T>C];[476T>C](;)1083A>C"
        ... ).description
        >>> desc.phase
        <AllelePhase.TRANS: 'trans'>
        >>> len(desc.unphased_alleles)
        1
        >>> len(desc.unphased_alleles[0].variants)
        1

        Variants *in cis* on a single protein allele:

        >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg;Asn594del]").description
        >>> len(tuple(desc))
        1
        >>> desc.allele_two is None
        True
        >>> desc.phase is None
        True
        >>> len(desc.allele_one.variants)
        2

        Two protein alleles *in trans*:

        >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg];[Ser68=]").description
        >>> desc.allele_two is not None
        True
        >>> desc.phase
        <AllelePhase.TRANS: 'trans'>
        >>> len(desc.allele_one.variants)
        1
        >>> len(desc.allele_two.variants)
        1
    """

    allele_one: Allele[VariantT]
    allele_two: Allele[VariantT] | None = None
    phase: AllelePhase | None = None
    alleles_unphased: tuple[Allele[VariantT], ...] = ()

    def __iter__(self) -> Iterator[Allele[VariantT]]:
        """Iterate over alleles in description order.

        Yields:
            (Allele[VariantT]): Alleles in description order: ``allele_one``,
                then ``allele_two`` when present, followed by any entries in
                ``alleles_unphased``.

        Notes:
            Iteration preserves the structural order of the HGVS allele-variant
            description. It does not infer or reorder alleles by phase.

        Examples:
            >>> from tinyhgvs import parse_hgvs
            >>> desc = parse_hgvs("NM_004006.2:c.[2376G>C];[3103del]").description
            >>> len(tuple(desc))
            2

            >>> desc = parse_hgvs(
            ...     "NM_004006.2:c.[296T>G;476T>C];[476T>C](;)1083A>C"
            ... ).description
            >>> len(tuple(desc))
            3

            >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg];[Ser68=]").description
            >>> len(tuple(desc))
            2

            >>> desc = parse_hgvs("p.[Ser68Arg];[Asn594del](;)0").description
            >>> len(tuple(desc))
            3
        """
        yield self.allele_one
        if self.allele_two is not None:
            yield self.allele_two
        yield from self.alleles_unphased

    @property
    def phased_alleles(
        self,
    ) -> tuple[Allele[VariantT], Allele[VariantT]] | None:
        """Return the established phased allele pair, if present.

        Returns:
            (tuple[Allele[VariantT], Allele[VariantT]] | None): The established
                phased allele pair as ``(allele_one, allele_two)`` when this
                description contains two established alleles with an explicit
                phase relationship; otherwise, ``None``.

        Notes:
            This property reports only the primary phased allele pair
            represented by ``allele_one`` and ``allele_two``. Alleles in
            ``alleles_unphased`` are not included.

        Examples:

            A single allele does not establish a phased pair:

            >>> from tinyhgvs import parse_hgvs
            >>> desc = parse_hgvs("NC_000001.11:g.[123G>A;345del]").description
            >>> desc.phased_alleles is None
            True

            Two alleles with established phase return a pair:

            >>> desc = parse_hgvs("NM_004006.2:c.[2376G>C];[2376=]").description
            >>> desc.phase
            <AllelePhase.TRANS: 'trans'>
            >>> pair = desc.phased_alleles
            >>> pair is not None
            True
            >>> len(pair[0].variants), len(pair[1].variants)
            (1, 1)

            Two alleles with uncertain phase do not return a pair:

            >>> desc = parse_hgvs("NC_000001.11:g.123G>A(;)345del").description
            >>> desc.phase
            <AllelePhase.UNCERTAIN: 'uncertain'>
            >>> pair = desc.phased_alleles
            >>> pair is None
            True

            Additional alleles with uncertain relation to the established pair:

            >>> desc = parse_hgvs(
            ...     "NM_004006.2:c.[296T>G;476T>C];[476T>C](;)1083A>C"
            ... ).description
            >>> pair = desc.phased_alleles
            >>> pair is not None
            True
            >>> len(desc.alleles_unphased)
            1

            Two protein alleles with known phase:

            >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg];[Ser68=]").description
            >>> desc.phase
            <AllelePhase.TRANS: 'trans'>
            >>> pair = desc.phased_alleles
            >>> pair is not None
            True
            >>> len(pair[0].variants), len(pair[1].variants)
            (1, 1)

            Two predicted protein alleles with unknown phase:

            >>> desc = parse_hgvs("NP_003997.1:p.(Ser73Arg)(;)(Asn103del)").description
            >>> desc.phase
            <AllelePhase.UNCERTAIN: 'uncertain'>
            >>> desc.phased_alleles is None
            True
            >>> len(desc.unphased_alleles)
            0
        """
        if self.phase is AllelePhase.TRANS and self.allele_two is not None:
            return (self.allele_one, self.allele_two)
        return None

    @property
    def unphased_alleles(self) -> tuple[Allele[VariantT], ...]:
        """Return alleles with uncertain relation to the allele state
        established by ``allele_one`` and ``allele_two``.

        Returns:
            (tuple[Allele[VariantT], ...]): Alleles written after the
                established allele state whose relation to that state is
                uncertain. Empty when not present.

        Notes:
            This property exposes the additional alleles stored in
            ``alleles_unphased``. It does not include ``allele_one`` or
            ``allele_two``.

        Examples:
            Uncertain phase between two primary alleles does not create an
            unphased tail:

            >>> from tinyhgvs import parse_hgvs
            >>> desc = parse_hgvs("NC_000001.11:g.123G>A(;)345del").description
            >>> len(desc.unphased_alleles)
            0

            One additional unphased allele to the established state:

            >>> desc = parse_hgvs(
            ...     "NM_004006.2:c.[296T>G];[476T>C](;)1083A>C"
            ... ).description
            >>> len(desc.unphased_alleles)
            1
            >>> len(desc.unphased_alleles[0].variants)
            1

            Multiple additions of unphased alleles to the established state:

            >>> desc = parse_hgvs(
            ...     "NM_004006.2:c.[296T>G];[476T>C](;)1083A>C(;)1406del"
            ... ).description
            >>> len(desc.unphased_alleles)
            2

            One additional unphased protein allele to the established state:

            >>> desc = parse_hgvs("p.[Ser68Arg];[Asn594del](;)0").description
            >>> len(desc.unphased_alleles)
            1
            >>> desc.unphased_alleles[0].variants[0].effect.kind
            'no_protein_produced'

            No additional unphased protein alleles:

            >>> desc = parse_hgvs("NP_003997.1:p.(Ser73Arg)(;)(Asn103del)").description
            >>> len(desc.unphased_alleles)
            0
        """
        return self.alleles_unphased


__all__ = [
    "Accession",
    "Allele",
    "AllelePhase",
    "AlleleVariant",
    "CoordinateSystem",
    "Interval",
    "PositionT",
    "ReferenceSpec",
    "VariantT",
]
