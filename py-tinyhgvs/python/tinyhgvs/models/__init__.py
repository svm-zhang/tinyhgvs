"""Public Python data model for parsed HGVS variants.

The package is split into:

- :mod:`tinyhgvs.models.shared` for shared reference and coordinate models
- :mod:`tinyhgvs.models.nucleotide` for nucleotide positions, edits, and variants
- :mod:`tinyhgvs.models.protein` for protein positions, effects, and variants

Type Aliases:
    VariantDescription: Tagged union for supported top-level variant payloads.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

from .nucleotide import (
    LiteralSequenceComponent,
    NucleotideDeletionEdit,
    NucleotideDeletionInsertionEdit,
    NucleotideDuplicationEdit,
    NucleotideEdit,
    NucleotideInsertionEdit,
    NucleotideInversionEdit,
    NucleotideNoChangeEdit,
    NucleotidePosition,
    NucleotidePositionAnchor,
    NucleotideSequence,
    NucleotideSequenceComponent,
    NucleotideSequenceSegment,
    NucleotideSubstitutionEdit,
    NucleotideVariant,
    RepeatSequenceComponent,
    SegmentSequenceComponent,
)
from .protein import (
    ProteinDeletionEdit,
    ProteinDeletionInsertionEdit,
    ProteinDuplicationEdit,
    ProteinEdit,
    ProteinEditEffect,
    ProteinEffect,
    ProteinInsertionEdit,
    ProteinNoChangeEdit,
    ProteinNoProteinProducedEffect,
    ProteinPosition,
    ProteinSequence,
    ProteinSubstitutionEdit,
    ProteinUnknownEdit,
    ProteinUnknownEffect,
    ProteinVariant,
)
from .shared import (
    CoordinateSystem,
    CurrentReferenceSource,
    OtherReferenceSource,
    Range,
    ReferenceSpec,
    SequenceId,
    SequenceKind,
    SequenceSource,
)

VariantDescription: TypeAlias = NucleotideVariant | ProteinVariant
"""Tagged union for supported top-level variant models:

- [`NucleotideVariant`][tinyhgvs.models.nucleotide.NucleotideVariant]
- [`ProteinVariant`][tinyhgvs.models.protein.ProteinVariant]
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
        >>> variant.reference.primary.raw
        'NM_004006.2'
        >>> variant.coordinate_system.value
        'c'
        >>> variant_description = variant.description
        >>> variant_location = variant_description.location
        >>> variant_location.start.position
        357
        >>> variant_location.start.offset
        1
        >>> variant_location.start.anchor
        <NucleotidePositionAnchor.COORDINATE: 'coordinate'>
        >>> variant_location.end is None
        True
        >>> variant_edit = variant_description.edit
        >>> variant_edit
        NucleotideSubstitutionEdit(reference='G', alternate='A', kind='substitution')
    """

    reference: ReferenceSpec | None
    coordinate_system: CoordinateSystem
    description: VariantDescription


__all__ = [
    "CoordinateSystem",
    "CurrentReferenceSource",
    "HgvsVariant",
    "LiteralSequenceComponent",
    "NucleotideDeletionEdit",
    "NucleotideDeletionInsertionEdit",
    "NucleotideDuplicationEdit",
    "NucleotideEdit",
    "NucleotideInsertionEdit",
    "NucleotideInversionEdit",
    "NucleotideNoChangeEdit",
    "NucleotidePosition",
    "NucleotidePositionAnchor",
    "NucleotideSequence",
    "NucleotideSequenceComponent",
    "NucleotideSequenceSegment",
    "NucleotideSubstitutionEdit",
    "NucleotideVariant",
    "OtherReferenceSource",
    "ProteinDeletionEdit",
    "ProteinDeletionInsertionEdit",
    "ProteinDuplicationEdit",
    "ProteinEdit",
    "ProteinEditEffect",
    "ProteinEffect",
    "ProteinInsertionEdit",
    "ProteinNoChangeEdit",
    "ProteinNoProteinProducedEffect",
    "ProteinPosition",
    "ProteinSequence",
    "ProteinSubstitutionEdit",
    "ProteinUnknownEdit",
    "ProteinUnknownEffect",
    "ProteinVariant",
    "Range",
    "ReferenceSpec",
    "RepeatSequenceComponent",
    "SegmentSequenceComponent",
    "SequenceId",
    "SequenceKind",
    "SequenceSource",
    "VariantDescription",
]
