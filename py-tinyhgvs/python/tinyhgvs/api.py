"""Public parsing entry points for tinyhgvs."""

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
        A fully typed `tinyhgvs.models.HgvsVariant` instance.

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
        >>> variant.description.reference
        'G'
        >>> variant.description.alternate
        'A'

        A 5' UTR substitution keeps its signed coordinate:

        >>> utr = parse_hgvs("NM_007373.4:c.-1C>T")
        >>> utr.description.location.start.coordinate
        -1
        >>> utr.description.location.start.is_five_prime_utr
        True

        An RNA repeat:

        >>> repeat = parse_hgvs("NM_004006.3:r.-124_-123[14]")
        >>> len(repeat.description.edit.sequence)
        1
        >>> repeat.description.edit.sequence[0].quantity.count
        14
        >>> repeat.description.edit.sequence[0].unit is None
        True

        Two alleles *in trans*:

        >>> allele = parse_hgvs("NM_004006.2:c.[2376G>C];[2376=]")
        >>> allele.description.phase
        <AllelePhase.TRANS: 'trans'>
        >>> len(tuple(allele.description.allele_one))
        1
        >>> len(tuple(allele.description.allele_two))
        1

        A predicted protein consequence:

        >>> protein = parse_hgvs("NP_003997.1:p.(Trp24Ter)")
        >>> protein.description.is_predicted
        True
        >>> protein.description.edit.location.start.residue
        'Trp'

        A predicted protein no-change outcome:

        >>> no_change = parse_hgvs("NP_003997.1:p.(Cys188=)")
        >>> no_change.description.is_no_change
        True
        >>> no_change.description.is_predicted
        True

        A protein frameshift variant (long-format):

        >>> frameshift = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
        >>> frameshift_edit = frameshift.description.edit
        >>> frameshift_edit.to_residue
        'Pro'
        >>> frameshift_edit.stop.ordinal
        23

        A protein insertion with unknown amino-acid content:

        >>> insertion = parse_hgvs("NP_003997.1:p.Arg78_Gly79insXaa[23]")
        >>> insertion.description.edit.count
        23

        Alternative protein consequences:

        >>> alternatives = parse_hgvs("NP_003997.1:p.(Gly23GlufsTer7^Gly23CysfsTer26)")
        >>> alternatives.description.has_alternatives
        True
        >>> len(alternatives.description.edits)
        2
    """
    from ._tinyhgvs import parse_hgvs as _parse_hgvs

    return _parse_hgvs(input)


__all__ = ["parse_hgvs"]
