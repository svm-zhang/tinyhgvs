"""Protein-focused Python model types for parsed HGVS variants.

Type Aliases:
    ProteinEdit: Tagged union for supported protein edit models.
    ProteinOutcome: Tagged union for supported protein outcome models.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TypeAlias

from .repeat import Repeat
from .core import Location, OutcomeCertainty


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
        >>> position = variant.description.edit.location.start
        >>> position.residue
        'Trp'
        >>> position.ordinal
        24
    """

    residue: str
    ordinal: int


_ResidueChange: TypeAlias = str | tuple[str, ...]


class ProteinExtensionTerminal(str, Enum):
    """Protein terminus toward which an extension variant extends.

    Attributes:
        N: Extension toward the N-terminus.
        C: Extension toward the C-terminus.

    Examples:
        N-terminal extension:
        >>> from tinyhgvs import ProteinExtensionTerminal, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.2:p.Met1ext-5")
        >>> variant.description.edit.to_terminal is ProteinExtensionTerminal.N
        True

        C-terminal extension:
        >>> variant = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")
        >>> variant.description.edit.to_terminal is ProteinExtensionTerminal.C
        True
    """

    N = "N"
    C = "C"


class ProteinFrameshiftStop:
    """Base class for protein frameshift stop-codon state."""

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
    """Short-form frameshift where stop codon information is omitted.

    Examples:
        >>> from tinyhgvs import OmittedProteinFrameshiftStop, parse_hgvs
        >>> variant = parse_hgvs("NP_0123456.1:p.Arg97fs")
        >>> stop = variant.description.edit.stop
        >>> isinstance(stop, OmittedProteinFrameshiftStop)
        True
    """

    @property
    def is_omitted(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class UnknownProteinFrameshiftStop(ProteinFrameshiftStop):
    """Frameshift where the stop codon is not encountered or not known.

    Examples:
        >>> from tinyhgvs import UnknownProteinFrameshiftStop, parse_hgvs
        >>> variant = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer?")
        >>> stop = variant.description.edit.stop
        >>> isinstance(stop, UnknownProteinFrameshiftStop)
        True
    """

    @property
    def is_unknown(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class KnownProteinFrameshiftStop(ProteinFrameshiftStop):
    """Frameshift with a known stop codon ordinal.

    Examples:
        >>> from tinyhgvs import KnownProteinFrameshiftStop, parse_hgvs
        >>> variant = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
        >>> stop = variant.description.edit.stop
        >>> isinstance(stop, KnownProteinFrameshiftStop)
        True
        >>> stop.ordinal
        23
    """

    ordinal: int

    @property
    def is_known(self) -> bool:
        return True


class ProteinEdit:
    """Base class for concrete protein edits."""

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
    """Protein substitution to another residue symbol.

    Examples:
        A tryptophan residue is replaced by a termination codon:
        >>> from tinyhgvs import ProteinSubstitution, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Trp24Ter")
        >>> edit = variant.description.edit
        >>> isinstance(edit, ProteinSubstitution)
        True
        >>> edit.to
        'Ter'

        Alternative consequences are represented by a tuple of residues:
        >>> variant = parse_hgvs("NP_003997.1:p.(Gly719Ala^Ser)")
        >>> edit = variant.description.edit
        >>> edit.to
        ('Ala', 'Ser')
        >>> edit.has_alternatives
        True
    """

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
    """Protein deletion.

    Examples:
        A deletion spanning residues Lys23 to Val25:
        >>> from tinyhgvs import ProteinDeletion, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.2:p.Lys23_Val25del")
        >>> edit = variant.description.edit
        >>> isinstance(edit, ProteinDeletion)
        True
        >>> edit.location.start.residue
        'Lys'
        >>> edit.location.end.residue
        'Val'
    """

    location: Location[ProteinCoordinate]

    @property
    def is_deletion(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinDuplication(ProteinEdit):
    """Protein duplication.

    Examples:
        >>> from tinyhgvs import ProteinDuplication, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Ser68_Arg70dup")
        >>> isinstance(variant.description.edit, ProteinDuplication)
        True
    """

    location: Location[ProteinCoordinate]

    @property
    def is_duplication(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinRepeat(ProteinEdit):
    """Protein repeat edit.

    Examples:
        >>> from tinyhgvs import ProteinRepeat, parse_hgvs
        >>> variant = parse_hgvs("NP_0123456.1:p.Arg65_Ser67[12]")
        >>> edit = variant.description.edit
        >>> isinstance(edit, ProteinRepeat)
        True
        >>> edit.repeat.quantity.count
        12
    """

    location: Location[ProteinCoordinate]
    repeat: Repeat

    @property
    def is_repeat(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinFrameshift(ProteinEdit):
    """Protein frameshift consequence.

    Examples:
        A short-form protein frameshift variant:
        >>> from tinyhgvs import OmittedProteinFrameshiftStop, ProteinFrameshift, parse_hgvs
        >>> short = parse_hgvs("NP_0123456.1:p.Arg97fs")
        >>> edit = short.description.edit
        >>> isinstance(edit, ProteinFrameshift)
        True
        >>> edit.to_residue is None
        True
        >>> isinstance(edit.stop, OmittedProteinFrameshiftStop)
        True

        A long-form protein frameshift variant:
        >>> long = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
        >>> edit = long.description.edit
        >>> edit.to_residue
        'Pro'
        >>> edit.stop.ordinal
        23

        A predicted long-form frameshift can have an unknown stop:
        >>> predicted = parse_hgvs("NP_0123456.1:p.(Arg97ProfsTer?)")
        >>> predicted.description.is_predicted
        True
        >>> predicted.description.edit.stop.is_unknown
        True
    """

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
    """Protein extension consequence.

    Examples:
        An N-terminal extension:
        >>> from tinyhgvs import ProteinExtensionTerminal, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.2:p.Met1ext-5")
        >>> edit = variant.description.edit
        >>> edit.to_terminal
        <ProteinExtensionTerminal.N: 'N'>
        >>> edit.to_residue is None
        True
        >>> edit.terminal_ordinal
        -5

        A C-terminal extension with known new stop:
        >>> variant = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")
        >>> edit = variant.description.edit
        >>> edit.to_terminal
        <ProteinExtensionTerminal.C: 'C'>
        >>> edit.to_residue
        'Gln'
        >>> edit.terminal_ordinal
        17

        A C-terminal extension with unknown new stop:
        >>> variant = parse_hgvs("NP_003997.2:p.Ter327ArgextTer?")
        >>> variant.description.edit.terminal_ordinal is None
        True
    """

    location: Location[ProteinCoordinate]
    to_terminal: ProteinExtensionTerminal
    to_residue: str | None
    terminal_ordinal: int | None

    @property
    def is_extension(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinDeletionInsertion(ProteinEdit):
    """Protein deletion-insertion.

    Examples:
        One residue is deleted and replaced by two amino acids:
        >>> from tinyhgvs import ProteinDeletionInsertion, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Cys28delinsTrpVal")
        >>> edit = variant.description.edit
        >>> isinstance(edit, ProteinDeletionInsertion)
        True
        >>> edit.sequence
        ('Trp', 'Val')
    """

    location: Location[ProteinCoordinate]
    sequence: tuple[str, ...]

    @property
    def is_deletion_insertion(self) -> bool:
        return True


class ProteinInsertion(ProteinEdit):
    """Base class for protein insertion edits."""

    pass


@dataclass(frozen=True, slots=True)
class KnownProteinInsertion(ProteinInsertion):
    """Protein insertion with a known amino-acid sequence.

    The public surface stores inserted residues directly as a tuple.

    Examples:
        Three amino acids are inserted between residues 2 and 3:
        >>> from tinyhgvs import KnownProteinInsertion, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Lys2_Gly3insGlnSerLys")
        >>> edit = variant.description.edit
        >>> isinstance(edit, KnownProteinInsertion)
        True
        >>> edit.sequence
        ('Gln', 'Ser', 'Lys')
    """

    location: Location[ProteinCoordinate]
    sequence: tuple[str, ...]

    @property
    def is_insertion(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class UnknownProteinInsertion(ProteinInsertion):
    """Protein insertion with unknown amino-acid content.

    Examples:
        A bare ``insXaa`` represents one unknown inserted amino acid:
        >>> from tinyhgvs import UnknownProteinInsertion, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Ser332_Ser333insXaa")
        >>> edit = variant.description.edit
        >>> isinstance(edit, UnknownProteinInsertion)
        True
        >>> edit.count
        1

        The unknown count can also be written explicitly:
        >>> variant = parse_hgvs("NP_003997.1:p.Arg78_Gly79insXaa[23]")
        >>> variant.description.edit.count
        23
    """

    location: Location[ProteinCoordinate]
    count: int

    @property
    def is_insertion(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class TerminatingProteinInsertion(ProteinInsertion):
    """Protein insertion containing a terminating residue.

    Examples:
        >>> from tinyhgvs import TerminatingProteinInsertion, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Gln746_Lys747ins*63")
        >>> edit = variant.description.edit
        >>> isinstance(edit, TerminatingProteinInsertion)
        True
        >>> edit.ordinal
        63
    """

    location: Location[ProteinCoordinate]
    ordinal: int

    @property
    def is_insertion(self) -> bool:
        return True


class ProteinOutcome:
    """Base class for protein outcomes.

    Protein outcomes cover produced edits, no-change outcomes, not-produced
    outcomes, and unknown outcomes.

    Examples:
        >>> from tinyhgvs import ProteinNoChange, ProteinProduced, parse_hgvs
        >>> isinstance(parse_hgvs("NP_003997.1:p.Trp24Ter").description, ProteinProduced)
        True
        >>> isinstance(parse_hgvs("NP_003997.1:p.Cys188=").description, ProteinNoChange)
        True
    """

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
    """Protein outcome where a concrete protein edit is produced.

    Examples:
        An observed protein consequence is not predicted:
        >>> from tinyhgvs import ProteinProduced, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Trp24Ter")
        >>> isinstance(variant.description, ProteinProduced)
        True
        >>> variant.description.is_predicted
        False

        A parenthesized protein consequence is predicted:
        >>> predicted = parse_hgvs("NP_003997.1:p.(Trp24Ter)")
        >>> predicted.description.is_predicted
        True
        >>> predicted.description.edit.location.start.residue
        'Trp'
    """

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
    """Protein outcome with alternative produced consequences.

    Examples:
        >>> from tinyhgvs import ProteinProducedAlternatives, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.(Gly23GlufsTer7^Gly23CysfsTer26)")
        >>> isinstance(variant.description, ProteinProducedAlternatives)
        True
        >>> len(variant.description.edits)
        2
        >>> variant.description.edits[0].is_frameshift
        True
    """

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
    """Protein outcome where no protein product is made.

    Examples:
        >>> from tinyhgvs import ProteinNotProduced, parse_hgvs
        >>> variant = parse_hgvs("LRG_199p1:p.0")
        >>> isinstance(variant.description, ProteinNotProduced)
        True
        >>> variant.description.is_not_produced
        True
    """

    certainty: OutcomeCertainty

    @property
    def is_not_produced(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED


@dataclass(frozen=True, slots=True)
class ProteinUnknown(ProteinOutcome):
    """Protein outcome ``p.?`` where the consequence is unknown.

    Examples:
        >>> from tinyhgvs import ProteinUnknown, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.?")
        >>> isinstance(variant.description, ProteinUnknown)
        True
    """

    @property
    def is_unknown(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class ProteinNoChange(ProteinOutcome):
    """Protein outcome with no amino-acid change.

    Examples:
        A site-specific no-change outcome:
        >>> from tinyhgvs import ProteinNoChange, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Cys188=")
        >>> isinstance(variant.description, ProteinNoChange)
        True
        >>> variant.description.location.start.residue
        'Cys'
        >>> variant.description.location.start.ordinal
        188
        >>> variant.description.is_no_change
        True

        A parenthesized site-specific no-change outcome is predicted:
        >>> predicted = parse_hgvs("NP_003997.1:p.(Cys188=)")
        >>> predicted.description.is_predicted
        True
        >>> predicted.description.location.start.ordinal
        188

        Whole-protein no-change has no site-specific location:
        >>> whole = parse_hgvs("NP_003997.1:p.(=)")
        >>> whole.description.location is None
        True
        >>> whole.description.is_predicted
        True
    """

    location: Location[ProteinCoordinate] | None
    certainty: OutcomeCertainty

    @property
    def is_no_change(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED


__all__ = [
    "KnownProteinFrameshiftStop",
    "KnownProteinInsertion",
    "OmittedProteinFrameshiftStop",
    "ProteinCoordinate",
    "ProteinEdit",
    "ProteinExtension",
    "ProteinExtensionTerminal",
    "ProteinFrameshift",
    "ProteinFrameshiftStop",
    "ProteinRepeat",
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
    "TerminatingProteinInsertion",
    "UnknownProteinFrameshiftStop",
    "UnknownProteinInsertion",
]
