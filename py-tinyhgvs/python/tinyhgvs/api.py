"""Public parsing entry points for :mod:`tinyhgvs`."""

from __future__ import annotations

from .models import HgvsVariant


def parse_hgvs(input: str) -> HgvsVariant:
    """Parse an HGVS string into the public Python data model.

    The parser trims leading and trailing whitespace, detects the HGVS
    coordinate type, and returns a typed model describing the reference,
    location, and edit.

    Args:
        input: HGVS expression to parse.

    Returns:
        A fully typed :class:`~tinyhgvs.models.HgvsVariant` instance.

    Raises:
        TinyHGVSError: If the input is invalid or belongs to a recognized but
            unsupported HGVS family.

    Examples:
        A splice-site coding DNA substitution:

        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> variant.coordinate_system.value
        'c'
        >>> variant.description.location.start.coordinate
        357
        >>> variant.description.location.start.offset
        1
        >>> variant.description.edit
        NucleotideSubstitutionEdit(reference='G', alternate='A', kind='substitution')

        A 5' UTR substitution keeps its signed coordinate:

        >>> utr = parse_hgvs("NM_007373.4:c.-1C>T")
        >>> utr.description.location.start.coordinate
        -1
        >>> utr.description.location.start.is_five_prime_utr
        True

        An exact RNA repeat:

        >>> repeat = parse_hgvs("NM_004006.3:r.-124_-123[14]")
        >>> len(repeat.description.edit.blocks)
        1
        >>> repeat.description.edit.blocks[0].count
        14
        >>> repeat.description.edit.blocks[0].unit is None
        True

        A predicted protein consequence:

        >>> protein = parse_hgvs("NP_003997.1:p.(Trp24Ter)")
        >>> protein.description.is_predicted
        True
        >>> protein.description.effect.location.start.residue
        'Trp'
    """
    from ._tinyhgvs import parse_hgvs as _parse_hgvs

    return _parse_hgvs(input)


__all__ = ["parse_hgvs"]
