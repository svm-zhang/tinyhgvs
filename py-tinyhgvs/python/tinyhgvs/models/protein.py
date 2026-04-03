"""Protein-focused Python model types for parsed HGVS variants.

Type Aliases:
    ProteinEdit: Tagged union for supported protein edit models.
    ProteinEffect: Tagged union for supported protein consequence models.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Literal, TypeAlias

from .shared import Interval


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


class ProteinSequenceOmittedEdit(str, Enum):
    """Protein edits whose altered amino-acid sequence is not written explicitly.

    Attributes:
        UNKNOWN: A protein change is expected, but the exact consequence is not known.
        NO_CHANGE: No protein change, written as ``=``.
        DELETION: Deletion of the stated amino-acid interval.
        DUPLICATION: Duplication of the stated amino-acid interval.

    Examples:
        A predicted but unspecified consequence at Met1:
        >>> from tinyhgvs import ProteinSequenceOmittedEdit, parse_hgvs
        >>> variant = parse_hgvs("LRG_199p1:p.(Met1?)")
        >>> variant.description.effect.edit is ProteinSequenceOmittedEdit.UNKNOWN
        True

        An explicitly unchanged residue:
        >>> variant = parse_hgvs("NP_003997.1:p.Cys188=")
        >>> variant.description.effect.edit
        <ProteinSequenceOmittedEdit.NO_CHANGE: 'no_change'>
    """

    UNKNOWN = "unknown"
    NO_CHANGE = "no_change"
    DELETION = "deletion"
    DUPLICATION = "duplication"


@dataclass(frozen=True, slots=True)
class ProteinSubstitutionEdit:
    """Protein substitution to another residue symbol.

    Attributes:
        to: Amino acid or stop symbol replacing the reference residue.
        kind: Edit kind.

    Examples:
        A tryptophan residue is replaced by a termination codon.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.Trp24Ter")
        >>> variant_edit = variant.description.effect.edit
        >>> variant_edit.to
        'Ter'
        >>> variant_edit.kind
        'substitution'
    """

    to: str
    kind: Literal["substitution"] = field(init=False, default="substitution")


@dataclass(frozen=True, slots=True)
class ProteinInsertionEdit:
    """Model describing protein insertion.

    Attributes:
        sequence: Inserted amino-acid sequence.
        kind: Edit kind.

    Examples:
        Three amino acids are inserted between residues 2 and 3.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("p.Lys2_Gly3insGlnSerLys")
        >>> variant_edit = variant.description.effect.edit
        >>> variant_edit.sequence.residues
        ('Gln', 'Ser', 'Lys')
        >>> variant_edit.kind
        'insertion'
    """

    sequence: ProteinSequence
    kind: Literal["insertion"] = field(init=False, default="insertion")


@dataclass(frozen=True, slots=True)
class ProteinDeletionInsertionEdit:
    """Model describing protein deletion-insertion.

    Attributes:
        sequence: Replacement amino-acid sequence.
        kind: Edit kind.

    Examples:
        One residue is deleted and replaced by two amino acids.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("p.Cys28delinsTrpVal")
        >>> variant_edit = variant.description.effect.edit
        >>> variant_edit.sequence.residues
        ('Trp', 'Val')
        >>> variant_edit.kind
        'deletion_insertion'
    """

    sequence: ProteinSequence
    kind: Literal["deletion_insertion"] = field(
        init=False,
        default="deletion_insertion",
    )


@dataclass(frozen=True, slots=True)
class ProteinRepeatEdit:
    """Model describing a top-level protein repeat variant.

    Attributes:
        count: Number of repeated amino-acid units.
        kind: Edit kind.

    Examples:
        A protein repeat variant with repeat unit coming from an interval:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NP_0123456.1:p.Arg65_Ser67[12]")
        >>> variant.description.effect.edit.count
        12
        >>> variant.description.effect.edit.kind
        'repeat'
    """

    count: int
    kind: Literal["repeat"] = field(init=False, default="repeat")


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

    N = "n"
    C = "c"


@dataclass(frozen=True, slots=True)
class ProteinExtensionEdit:
    """Model describing a protein extension consequence.

    Attributes:
        to_terminal: Protein terminus toward which an extension variant extends.
        to_residue: Residue replacing the reference stop codon in a C-terminal
            extension. None for N-terminal extension.
        terminal_ordinal: New terminal ordinal. Negative for N-terminal
            extension, positive for C-terminal extension with known stop, and
            None when the new stop is unknown.
        kind: Edit kind.

    Examples:
        An N-terminal extension:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NP_003997.2:p.Met1ext-5")
        >>> variant_edit = variant.description.effect.edit
        >>> variant_edit.to_terminal
        <ProteinExtensionTerminal.N: 'n'>
        >>> variant_edit.to_residue is None
        True
        >>> variant_edit.terminal_ordinal
        -5

        A C-terminal extension with known new stop:
        >>> variant = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")
        >>> variant_edit = variant.description.effect.edit
        >>> variant_edit.to_terminal
        <ProteinExtensionTerminal.C: 'c'>
        >>> variant_edit.to_residue
        'Gln'
        >>> variant_edit.terminal_ordinal
        17

        A C-terminal extension with unknown new stop:
        >>> variant = parse_hgvs("p.Ter327ArgextTer?")
        >>> variant.description.effect.edit.terminal_ordinal is None
        True
    """

    to_terminal: ProteinExtensionTerminal
    to_residue: str | None
    terminal_ordinal: int | None
    kind: Literal["extension"] = field(init=False, default="extension")


class ProteinFrameshiftStopKind(str, Enum):
    """Model describing a stop codon is known (long-form), or omitted (short-form),
    or unknown (not encountered) due to a frameshift event.

    Attributes:
        OMITTED: Short-form frameshift where stop codon information is omitted.
        UNKNOWN: Long-form frameshift with unknown stop codon, such as ``p.Arg97ProfsTer?``.
        KNOWN: Long-form frameshift with known stop codon, such as ``p.Arg97ProfsTer23``.

    Examples:
        Protein frameshift variant with unknown stop codon:
        >>> from tinyhgvs import parse_hgvs
        >>> unknown = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer?")
        >>> unknown.description.effect.edit.stop.kind
        <ProteinFrameshiftStopKind.UNKNOWN: 'unknown'>

        Short-form protein frameshift variant:
        >>> short = parse_hgvs("NP_0123456.1:p.Arg97fs")
        >>> short.description.effect.edit.stop.kind
        <ProteinFrameshiftStopKind.OMITTED: 'omitted'>

    """

    OMITTED = "omitted"
    UNKNOWN = "unknown"
    KNOWN = "known"


@dataclass(frozen=True, slots=True)
class ProteinFrameshiftStop:
    """Model describing stop codon information in a protein frameshift edit.

    Attributes:
        ordinal: Stop codon ordinal. None when stop codon information is either
            omitted (short-form) or unknown.
        kind: Whether the stop is omitted, unknown, or known.

    Examples:
        Protein frameshift variant with a known stop codon:
        >>> from tinyhgvs import parse_hgvs
        >>> long = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
        >>> long_stop = long.description.effect.edit.stop
        >>> long_stop.ordinal
        23
        >>> long_stop.kind
        <ProteinFrameshiftStopKind.KNOWN: 'known'>

        Short-form frameshift variant leaves the stop ordinal optional:
        >>> short = parse_hgvs("NP_0123456.1:p.Arg97fs")
        >>> short_stop = short.description.effect.edit.stop
        >>> short_stop.ordinal is None
        True
        >>> short_stop.kind
        <ProteinFrameshiftStopKind.OMITTED: 'omitted'>

    """

    ordinal: int | None
    kind: ProteinFrameshiftStopKind


@dataclass(frozen=True, slots=True)
class ProteinFrameshiftEdit:
    """Model describing a protein frameshift consequence.

    Attributes:
        to_residue: First newly encoded residue. None when in short-form.
        stop: Model describing stop codon information in frameshift variant.
        kind: Edit kind.

    Examples:
        A short-form protein frameshift variant:
        >>> from tinyhgvs import parse_hgvs
        >>> short = parse_hgvs("NP_0123456.1:p.Arg97fs")
        >>> short_edit = short.description.effect.edit
        >>> short_edit.to_residue is None
        True
        >>> short_edit.stop.kind
        <ProteinFrameshiftStopKind.OMITTED: 'omitted'>

        A long-form protein frameshift variants:
        >>> long = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
        >>> long_edit = long.description.effect.edit
        >>> long_edit.to_residue
        'Pro'
        >>> long_edit.stop.ordinal
        23

        A predicted long-form protein frameshift variant:
        >>> predicted = parse_hgvs("p.(Arg97ProfsTer?)")
        >>> predicted.description.is_predicted
        True
        >>> predicted.description.effect.edit.stop.kind
        <ProteinFrameshiftStopKind.UNKNOWN: 'unknown'>
    """

    to_residue: str | None
    stop: ProteinFrameshiftStop
    kind: Literal["frameshift"] = field(init=False, default="frameshift")


ProteinEdit: TypeAlias = (
    ProteinSequenceOmittedEdit
    | ProteinSubstitutionEdit
    | ProteinInsertionEdit
    | ProteinDeletionInsertionEdit
    | ProteinRepeatEdit
    | ProteinExtensionEdit
    | ProteinFrameshiftEdit
)
"""Tagged union for supported protein edit models:

- [`ProteinSequenceOmittedEdit`][tinyhgvs.models.protein.ProteinSequenceOmittedEdit]
- [`ProteinSubstitutionEdit`][tinyhgvs.models.protein.ProteinSubstitutionEdit]
- [`ProteinInsertionEdit`][tinyhgvs.models.protein.ProteinInsertionEdit]
- [`ProteinDeletionInsertionEdit`][tinyhgvs.models.protein.ProteinDeletionInsertionEdit]
- [`ProteinRepeatEdit`][tinyhgvs.models.protein.ProteinRepeatEdit]
- [`ProteinExtensionEdit`][tinyhgvs.models.protein.ProteinExtensionEdit]
- [`ProteinFrameshiftEdit`][tinyhgvs.models.protein.ProteinFrameshiftEdit]
"""


@dataclass(frozen=True, slots=True)
class ProteinUnknownEffect:
    """Model describing the protein consequence ``p.?``.

    Attributes:
        kind: Effect kind.

    Examples:
        The protein consequence is entirely unknown.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NP_003997.1:p.?")
        >>> variant.description.effect.kind
        'unknown'
    """

    kind: Literal["unknown"] = field(init=False, default="unknown")


@dataclass(frozen=True, slots=True)
class ProteinNoProteinProducedEffect:
    """Model describing the protein consequence ``p.0``.

    Attributes:
        kind: Effect kind.

    Examples:
        The variant predicts that no protein product is made.
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("LRG_199p1:p.0")
        >>> variant.description.effect.kind
        'no_protein_produced'
    """

    kind: Literal["no_protein_produced"] = field(
        init=False,
        default="no_protein_produced",
    )


@dataclass(frozen=True, slots=True)
class ProteinEditEffect:
    """Concrete protein consequence at a known amino-acid interval.

    Attributes:
        location: Amino-acid range where the edit occurs.
        edit: Protein edit applied at that range.
        kind: Effect kind.

    Examples:
        A deletion spanning residues Lys23 to Val25 is represented by a protein
        location and a protein deletion edit.
        >>> from tinyhgvs import ProteinSequenceOmittedEdit, parse_hgvs
        >>> variant = parse_hgvs("NP_003997.2:p.Lys23_Val25del")
        >>> effect = variant.description.effect
        >>> effect.location.start.residue
        'Lys'
        >>> effect.location.end.residue
        'Val'
        >>> effect.edit is ProteinSequenceOmittedEdit.DELETION
        True

        A protein frameshift consequence:
        >>> variant = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
        >>> effect = variant.description.effect
        >>> effect.location.start.residue
        'Arg'
        >>> effect.edit.to_residue
        'Pro'

        A protein extension consequence:
        >>> variant = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")
        >>> effect = variant.description.effect
        >>> effect.location.start.residue
        'Ter'
        >>> effect.edit.kind
        'extension'
        >>> effect.edit.terminal_ordinal
        17
    """

    location: Interval[ProteinCoordinate]
    edit: ProteinEdit
    kind: Literal["edit"] = field(init=False, default="edit")


ProteinEffect: TypeAlias = (
    ProteinUnknownEffect | ProteinNoProteinProducedEffect | ProteinEditEffect
)
"""Tagged union for supported protein consequence models:

- [`ProteinUnknownEffect`][tinyhgvs.models.protein.ProteinUnknownEffect]
- [`ProteinNoProteinProducedEffect`][tinyhgvs.models.protein.ProteinNoProteinProducedEffect]
- [`ProteinEditEffect`][tinyhgvs.models.protein.ProteinEditEffect]
"""


@dataclass(frozen=True, slots=True)
class ProteinVariant:
    """Parsed protein-level consequence.

    Attributes:
        is_predicted: Whether the effect was written in parentheses.
        effect: Parsed protein consequence model.

    Examples:
        An observed protein consequence is not predicted:
        >>> from tinyhgvs import parse_hgvs
        >>> parse_hgvs("NP_003997.1:p.Trp24Ter").description.is_predicted
        False

        A parenthesized protein consequence is predicted:
        >>> parse_hgvs("LRG_199p1:p.(Met1?)").description.is_predicted
        True

        Frameshift consequences:
        >>> variant = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
        >>> variant.description.effect.edit.kind
        'frameshift'

        Extension consequences:
        >>> variant = parse_hgvs("NP_003997.2:p.Ter110GlnextTer17")
        >>> variant.description.effect.edit.kind
        'extension'
    """

    is_predicted: bool
    effect: ProteinEffect


ProteinPosition = ProteinCoordinate


__all__ = [
    "ProteinCoordinate",
    "ProteinDeletionInsertionEdit",
    "ProteinEdit",
    "ProteinEditEffect",
    "ProteinEffect",
    "ProteinExtensionEdit",
    "ProteinExtensionTerminal",
    "ProteinFrameshiftEdit",
    "ProteinFrameshiftStop",
    "ProteinFrameshiftStopKind",
    "ProteinInsertionEdit",
    "ProteinNoProteinProducedEffect",
    "ProteinRepeatEdit",
    "ProteinSequence",
    "ProteinSequenceOmittedEdit",
    "ProteinSubstitutionEdit",
    "ProteinUnknownEffect",
    "ProteinVariant",
]
