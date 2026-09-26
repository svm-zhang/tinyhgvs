from dataclasses import dataclass
from enum import Enum
from typing import Generic, Iterator, TypeVar

_T = TypeVar("_T")


class AlleleStateCertainty(str, Enum):
    CERTAIN = "certain"
    UNCERTAIN = "uncertain"


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
class Allele(Generic[_T]):
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

    variants: tuple[_T, ...]
    state_certainty: AlleleStateCertainty

    @property
    def is_uncertain(self) -> bool:
        return self.state_certainty is AlleleStateCertainty.UNCERTAIN

    def __iter__(self) -> Iterator[_T]:
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


class AlleleForm(Generic[_T]):
    """Base class for top-level allele forms.

    Allele forms describe how allele-level syntax is written at the outer
    level. A single allele form covers ordinary allele syntax, while derived
    and alternative forms cover comma and ``^`` syntax.

    Examples:
        Ordinary allele syntax:
        >>> from tinyhgvs import AlleleVariant, parse_hgvs
        >>> desc = parse_hgvs("NM_004006.2:c.[2376G>C];[2376=]").description
        >>> isinstance(desc, AlleleVariant)
        True
        >>> desc.is_single
        True

        Derived allele syntax:
        >>> from tinyhgvs import DerivedAlleleForm
        >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg,Asn594del]").description
        >>> isinstance(desc, DerivedAlleleForm)
        True
        >>> desc.is_derived
        True

        Alternative allele syntax:
        >>> from tinyhgvs import AlternativeAlleleForm
        >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg]^[Asn594del]").description
        >>> isinstance(desc, AlternativeAlleleForm)
        True
        >>> desc.is_alternative
        True
    """

    @property
    def is_single(self) -> bool:
        return False

    @property
    def is_derived(self) -> bool:
        return False

    @property
    def is_alternative(self) -> bool:
        return False


@dataclass(frozen=True, slots=True)
class AlleleVariant(AlleleForm[_T]):
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
        unphased: Later variants written in uncertain relation to the
            established allele state.

    Examples:
        Variants *in cis* on a single allele:

        >>> from tinyhgvs import parse_hgvs
        >>> desc = parse_hgvs("NC_000023.10:g.[30683643A>G;33038273T>G]").description
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
        ...     "NC_000001.11:g.[123G>A];[345del](;)789dup"
        ... ).description
        >>> desc.phase
        <AllelePhase.TRANS: 'trans'>
        >>> len(desc.unphased)
        1
        >>> desc.unphased[0].location.start.coordinate
        789

        Variants *in cis* on a single protein allele:

        >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg;Asn594del]").description
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

    allele_one: Allele[_T]
    allele_two: Allele[_T] | None
    phase: AllelePhase | None
    unphased: tuple[_T, ...]

    @property
    def is_single(self) -> bool:
        return True

    @property
    def has_second_allele(self) -> bool:
        return self.allele_two is not None

    @property
    def has_unphased_variants(self) -> bool:
        return bool(self.unphased)

    @property
    def phased_alleles(
        self,
    ) -> tuple[Allele[_T], Allele[_T]] | None:
        """Return the established phased allele pair, if present.

        Returns:
            (tuple[Allele[VariantT], Allele[VariantT]] | None): The established
                phased allele pair as ``(allele_one, allele_two)`` when this
                description contains two established alleles with an explicit
                phase relationship; otherwise, ``None``.

        Notes:
            This property reports only the primary phased allele pair
            represented by ``allele_one`` and ``allele_two``. Alleles in
            ``unphased`` are not included.

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
            ...     "NC_000001.11:g.[123G>A];[345del](;)789dup"
            ... ).description
            >>> pair = desc.phased_alleles
            >>> pair is not None
            True
            >>> len(desc.unphased)
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
            >>> len(desc.unphased)
            0
        """
        if self.phase is AllelePhase.TRANS and self.allele_two is not None:
            return (self.allele_one, self.allele_two)
        return None


@dataclass(frozen=True, slots=True)
class DerivedAlleleForm(AlleleForm[_T]):
    """Allele form for derived outcomes written with comma syntax.

    Examples:
        >>> from tinyhgvs import DerivedAlleleForm, parse_hgvs
        >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg,Asn594del]").description
        >>> isinstance(desc, DerivedAlleleForm)
        True
        >>> len(desc.outcomes)
        2
        >>> desc.outcomes[0].edit.to
        'Arg'
    """

    outcomes: tuple[_T, ...]

    @property
    def is_derived(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class AlternativeAlleleForm(AlleleForm[_T]):
    """Allele form for alternative allele variants written with ``^`` syntax.

    Examples:
        >>> from tinyhgvs import AlternativeAlleleForm, parse_hgvs
        >>> desc = parse_hgvs("NP_003997.1:p.[Ser68Arg]^[Asn594del]").description
        >>> isinstance(desc, AlternativeAlleleForm)
        True
        >>> len(desc.alternatives)
        2
        >>> desc.alternatives[0].allele_one.variants[0].edit.to
        'Arg'
    """

    alternatives: tuple[AlleleVariant[_T], ...]

    @property
    def is_alternative(self) -> bool:
        return True


__all__ = [
    "Allele",
    "AlleleStateCertainty",
    "AllelePhase",
    "AlleleVariant",
    "AlleleForm",
    "DerivedAlleleForm",
    "AlternativeAlleleForm",
]
