from dataclasses import dataclass
from typing import TypeAlias

from .nucleotide import NucleotideEdit


@dataclass(frozen=True, slots=True)
class CodingDnaUnknown:
    """Coding-DNA outcome ``c.?`` where the consequence is unknown.

    Examples:
        >>> from tinyhgvs import CodingDnaUnknown, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.?")
        >>> isinstance(variant.description, CodingDnaUnknown)
        True
    """

    pass


CodingDnaOutcome: TypeAlias = NucleotideEdit | CodingDnaUnknown
"""Tagged union for supported coding-DNA outcomes:

- [`NucleotideEdit`][tinyhgvs.models.nucleotide.NucleotideEdit]
- [`CodingDnaUnknown`][tinyhgvs.models.cdna.CodingDnaUnknown]
"""

__all__ = [
    "CodingDnaUnknown",
    "CodingDnaOutcome",
]
