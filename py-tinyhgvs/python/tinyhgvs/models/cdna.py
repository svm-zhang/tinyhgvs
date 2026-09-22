from dataclasses import dataclass
from typing import TypeAlias

from .nucleotide import NucleotideEdit


@dataclass(frozen=True, slots=True)
class CodingDnaUnknown:
    pass


CodingDnaOutcome: TypeAlias = NucleotideEdit | CodingDnaUnknown

__all__ = [
    "CodingDnaUnknown",
    "CodingDnaOutcome",
]
