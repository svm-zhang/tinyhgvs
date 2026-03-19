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

ProteinEdit: TypeAlias = (
    ProteinSequenceOmittedEdit
    | ProteinSubstitutionEdit
    | ProteinInsertionEdit
    | ProteinDeletionInsertionEdit
)
"""Tagged union for supported protein edit models:

- [`ProteinSequenceOmittedEdit`][tinyhgvs.models.protein.ProteinSequenceOmittedEdit]
- [`ProteinSubstitutionEdit`][tinyhgvs.models.protein.ProteinSubstitutionEdit]
- [`ProteinInsertionEdit`][tinyhgvs.models.protein.ProteinInsertionEdit]
- [`ProteinDeletionInsertionEdit`][tinyhgvs.models.protein.ProteinDeletionInsertionEdit]
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
    "ProteinInsertionEdit",
    "ProteinNoProteinProducedEffect",
    "ProteinSequence",
    "ProteinSequenceOmittedEdit",
    "ProteinSubstitutionEdit",
    "ProteinUnknownEffect",
    "ProteinVariant",
]
