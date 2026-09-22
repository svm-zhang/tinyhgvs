"""Protein-focused Python model types for parsed HGVS variants.

Type Aliases:
    ProteinEdit: Tagged union for supported protein edit models.
    ProteinEffect: Tagged union for supported protein consequence models.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TypeAlias

from .shared import Location, OutcomeCertainty, Repeat


@dataclass(frozen=True, slots=True)
class ProteinCoordinate:
    """Protein position written as residue symbol plus ordinal.

    Attributes:
        residue: Amino-acid symbol.
        ordinal: Amino-acid position.

    Examples:
        A protein substitution at residue 24 is located using the amino-acid
        symbol and ordinal together.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Trp24Ter")
        >>> position = variant.description.effect.location.start
        >>> position.residue
        'Trp'
        >>> position.ordinal
        24
    """

    residue: str
    ordinal: int


_ResidueChange: TypeAlias = str | tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ProteinSequence:
    """Ordered amino-acid sequence used by insertions and deletion-insertions.

    Attributes:
        residues: Ordered tuple of amino-acid symbols.

    Examples:
        A protein insertion adds three amino acids in order.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("p.Lys2_Gly3insGlnSerLys")
        >>> variant_edit = variant.description.effect.edit
        >>> variant_edit.sequence.residues
        ('Gln', 'Ser', 'Lys')
    """

    residues: tuple[str, ...]


class ProteinExtensionTerminal(str, Enum):
    """Protein terminus toward which an extension variant extends.

    Attributes:
        N: Extension toward the N-terminus.
        C: Extension toward the C-terminus.

    Examples:
        N-terminal extension:
        >>> from tinyhgvs import ProteinExtensionTerminal, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.2:p.Met1ext-5")
        >>> variant.description.effect.edit.to_terminal is ProteinExtensionTerminal.N
        True

        C-terminal extension:
        >>> variant = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")
        >>> variant.description.effect.edit.to_terminal is ProteinExtensionTerminal.C
        True
    """

    N = "N"
    C = "C"


class ProteinFrameshiftStop:
    @property
    def is_omitted(self) -> bool:
        return False

    @property
    def is_unknown(self) -> bool:
        return False

    @property
    def is_known(self) -> bool:
        return False


@dataclass(frozen=True, slots=True)
class OmittedProteinFrameshiftStop(ProteinFrameshiftStop):
    @property
    def is_omitted(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class UnknownProteinFrameshiftStop(ProteinFrameshiftStop):
    @property
    def is_unknown(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class KnownProteinFrameshiftStop(ProteinFrameshiftStop):
    ordinal: int

    @property
    def is_known(self) -> bool:
        return True


class ProteinEdit:
    @property
    def is_substitution(self) -> bool:
        return False

    @property
    def is_deletion(self) -> bool:
        return False

    @property
    def is_duplication(self) -> bool:
        return False

    @property
    def is_repeat(self) -> bool:
        return False

    @property
    def is_extension(self) -> bool:
        return False

    @property
    def is_frameshift(self) -> bool:
        return False

    @property
    def is_insertion(self) -> bool:
        return False

    @property
    def is_deletion_insertion(self) -> bool:
        return False


@dataclass(frozen=True, slots=True)
class ProteinSubstitution(ProteinEdit):
    location: Location[ProteinCoordinate]
    to: _ResidueChange

    @property
    def is_substitution(self) -> bool:
        return True

    @property
    def has_alternatives(self) -> bool:
        return isinstance(self.to, tuple)


@dataclass(frozen=True, slots=True)
class ProteinDeletion(ProteinEdit):
    location: Location[ProteinCoordinate]

    @property
    def is_deletion(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinDuplication(ProteinEdit):
    location: Location[ProteinCoordinate]

    @property
    def is_duplication(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinRepeat(ProteinEdit):
    location: Location[ProteinCoordinate]
    repeat: Repeat

    @property
    def is_repeat(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinFrameshift(ProteinEdit):
    location: Location[ProteinCoordinate]
    to_residue: _ResidueChange | None
    stop: ProteinFrameshiftStop

    @property
    def is_frameshift(self) -> bool:
        return True

    @property
    def has_alternative_residues(self) -> bool:
        return isinstance(self.to_residue, tuple)


@dataclass(frozen=True, slots=True)
class ProteinExtension(ProteinEdit):
    location: Location[ProteinCoordinate]
    to_terminal: ProteinExtensionTerminal
    to_residue: str | None
    terminal_ordinal: int | None

    @property
    def is_extension(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinDeletionInsertion(ProteinEdit):
    location: Location[ProteinCoordinate]
    sequence: tuple[str, ...]

    @property
    def is_deletion_insertion(self) -> bool:
        return True


class ProteinInsertion(ProteinEdit):
    pass


@dataclass(frozen=True, slots=True)
class KnownProteinInsertion(ProteinInsertion):
    location: Location[ProteinCoordinate]
    sequence: tuple[str, ...]

    @property
    def is_insertion(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class UnknownProteinInsertion(ProteinInsertion):
    location: Location[ProteinCoordinate]
    count: int

    @property
    def is_insertion(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class TerminatingProteinInsertion(ProteinInsertion):
    location: Location[ProteinCoordinate]
    ordinal: int

    @property
    def is_insertion(self) -> bool:
        return True


class ProteinOutcome:
    @property
    def is_produced(self) -> bool:
        return False

    @property
    def is_no_change(self) -> bool:
        return False

    @property
    def is_not_produced(self) -> bool:
        return False

    @property
    def is_unknown(self) -> bool:
        return False


@dataclass(frozen=True, slots=True)
class ProteinProduced(ProteinOutcome):
    edit: ProteinEdit
    certainty: OutcomeCertainty

    @property
    def is_produced(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED

    @property
    def has_alternatives(self) -> bool:
        return False


@dataclass(frozen=True, slots=True)
class ProteinProducedAlternatives(ProteinOutcome):
    edits: tuple[ProteinEdit, ...]
    certainty: OutcomeCertainty

    @property
    def is_produced(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED

    @property
    def has_alternatives(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinNotProduced(ProteinOutcome):
    certainty: OutcomeCertainty

    @property
    def is_not_produced(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED


@dataclass(frozen=True, slots=True)
class ProteinUnknown(ProteinOutcome):
    @property
    def is_unknown(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinNoChange(ProteinOutcome):
    certainty: OutcomeCertainty

    @property
    def is_no_change(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED


__all__ = [
    "ProteinCoordinate",
    "ProteinEdit",
    "ProteinExtension",
    "ProteinExtensionTerminal",
    "ProteinFrameshift",
    "ProteinFrameshiftStop",
    "ProteinRepeat",
    "ProteinSequence",
    "ProteinSubstitution",
    "ProteinDeletion",
    "ProteinInsertion",
    "ProteinDeletionInsertion",
    "ProteinDuplication",
    "ProteinOutcome",
    "ProteinNoChange",
    "ProteinProduced",
    "ProteinProducedAlternatives",
    "ProteinNotProduced",
    "ProteinUnknown",
]
