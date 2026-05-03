"""Public Python data model for parsed HGVS variants.

The package is split into:

- :mod:`tinyhgvs.models.shared` for shared reference and coordinate models
- :mod:`tinyhgvs.models.nucleotide` for nucleotide coordinates, edits, and variants
- :mod:`tinyhgvs.models.protein` for protein coordinates, effects, and variants

Type Aliases:
    VariantDescription: Tagged union for supported top-level variant models.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

from .nucleotide import (
    CopiedSequenceItem,
    LiteralSequenceItem,
    NucleotideAnchor,
    NucleotideCoordinate,
    NucleotideDeletionInsertionEdit,
    NucleotideEdit,
    NucleotideInsertionEdit,
    NucleotideRepeatBlock,
    NucleotideRepeatEdit,
    NucleotideSequenceItem,
    NucleotideSequenceOmittedEdit,
    NucleotideSubstitutionEdit,
    NucleotideVariant,
    RepeatSequenceItem,
)
from .protein import (
    ProteinCoordinate,
    ProteinDeletionInsertionEdit,
    ProteinEdit,
    ProteinEditEffect,
    ProteinEffect,
    ProteinExtensionEdit,
    ProteinExtensionTerminal,
    ProteinFrameshiftEdit,
    ProteinFrameshiftStop,
    ProteinFrameshiftStopKind,
    ProteinInsertionEdit,
    ProteinNoProteinProducedEffect,
    ProteinRepeatEdit,
    ProteinSequence,
    ProteinSequenceOmittedEdit,
    ProteinSubstitutionEdit,
    ProteinUnknownEffect,
    ProteinVariant,
)
from .shared import (
    Accession,
    Allele,
    AllelePhase,
    AlleleVariant,
    CoordinateSystem,
    Interval,
    Location,
    ReferenceSpec,
)

VariantDescription: TypeAlias = (
    NucleotideVariant
    | AlleleVariant[NucleotideVariant]
    | ProteinVariant
    | AlleleVariant[ProteinVariant]
)
"""Tagged union for supported top-level variant models:

- [`NucleotideVariant`][tinyhgvs.models.nucleotide.NucleotideVariant]
- `AlleleVariant[NucleotideVariant]`
- [`ProteinVariant`][tinyhgvs.models.protein.ProteinVariant]
- `AlleleVariant[ProteinVariant]`
"""


@dataclass(frozen=True, slots=True)
class HgvsVariant:
    """Top-level model describing a parsed HGVS variant.

    This is the root object returned by :func:`tinyhgvs.parse_hgvs`. It ties the
    reference field, coordinate system, and parsed variant description together.

    Attributes:
        reference: Reference sequence field preceding the ``:`` when present.
        coordinate_system: HGVS coordinate type.
        description: Model describing a nucleotide or protein variant.

    Examples:
        A coding DNA splice-site substitution parsed into reference, location,
        and edit models:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NM_004006.2:c.357+1G>A")
        >>> variant.reference.primary.id
        'NM_004006.2'
        >>> variant.coordinate_system.value
        'c'
        >>> variant_description = variant.description
        >>> variant_location = variant_description.location
        >>> variant_location.start.coordinate
        357
        >>> variant_location.start.offset
        1
        >>> variant_location.start.anchor
        <NucleotideAnchor.ABSOLUTE: 'absolute'>
        >>> variant_location.end is None
        True
        >>> variant_edit = variant_description.edit
        >>> variant_edit
        NucleotideSubstitutionEdit(reference='G', alternate='A', kind='substitution')

        A 5' UTR substitution keeps the signed coordinate from the HGVS string:
        >>> utr = parse_hgvs("NM_007373.4:c.-1C>T")
        >>> utr.description.location.start.coordinate
        -1
        >>> utr.description.location.start.is_five_prime_utr
        True

        A protein frameshift is still exposed through the same top-level
        variant container:
        >>> protein = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
        >>> protein.description.effect.edit.kind
        'frameshift'
    """

    reference: ReferenceSpec | None
    coordinate_system: CoordinateSystem
    description: VariantDescription


__all__ = [
    "Accession",
    "Allele",
    "AllelePhase",
    "AlleleVariant",
    "CopiedSequenceItem",
    "CoordinateSystem",
    "HgvsVariant",
    "Interval",
    "Location",
    "LiteralSequenceItem",
    "NucleotideDeletionInsertionEdit",
    "NucleotideAnchor",
    "NucleotideCoordinate",
    "NucleotideEdit",
    "NucleotideInsertionEdit",
    "NucleotideRepeatBlock",
    "NucleotideRepeatEdit",
    "NucleotideSequenceItem",
    "NucleotideSequenceOmittedEdit",
    "NucleotideSubstitutionEdit",
    "NucleotideVariant",
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
    "ReferenceSpec",
    "RepeatSequenceItem",
    "VariantDescription",
]
