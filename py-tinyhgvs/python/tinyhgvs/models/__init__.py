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

from .allele import (
    Allele,
    AlleleForm,
    AllelePhase,
    AlleleStateCertainty,
    AlleleVariant,
    AlternativeAlleleForm,
    DerivedAlleleForm,
)
from .cdna import (
    CodingDnaOutcome,
    CodingDnaUnknown,
)
from .nucleotide import (
    CopiedSequence,
    LiteralSequence,
    NucleotideAnchor,
    NucleotideCoordinate,
    NucleotideDeletion,
    NucleotideDeletionInsertion,
    NucleotideDuplication,
    NucleotideEdit,
    NucleotideInsertion,
    NucleotideInversion,
    NucleotideNoChange,
    NucleotideRepeat,
    NucleotideRepeatBlock,
    NucleotideSequenceItem,
    NucleotideSubstitution,
    RepeatSequenceItem,
)
from .protein import (
    KnownProteinFrameshiftStop,
    OmittedProteinFrameshiftStop,
    ProteinCoordinate,
    ProteinDeletion,
    ProteinDeletionInsertion,
    ProteinDuplication,
    ProteinEdit,
    ProteinExtension,
    ProteinExtensionTerminal,
    ProteinFrameshift,
    ProteinFrameshiftStop,
    ProteinInsertion,
    ProteinNoChange,
    ProteinNotProduced,
    ProteinOutcome,
    ProteinProduced,
    ProteinProducedAlternatives,
    ProteinRepeat,
    ProteinSequence,
    ProteinSubstitution,
    ProteinUnknown,
    UnknownProteinFrameshiftStop,
)
from .repeat import (
    KnownQuantity,
    KnownRepeatUnit,
    Repeat,
    RepeatUnit,
    UncertainQuantity,
    UnknownQuantity,
    UnknownRepeatUnit,
)
from .rna import (
    RnaIndeterminate,
    RnaNoChange,
    RnaNotProduced,
    RnaOutcome,
    RnaProduced,
    RnaUncertainSplicing,
    RnaUnknown,
)
from .shared import (
    Accession,
    CoordinateSystem,
    KnownLocation,
    Location,
    OutcomeCertainty,
    PossibleRange,
    Quantity,
    ReferenceSpec,
    UncertainLocation,
)

_VariantDescription: TypeAlias = (
    NucleotideEdit
    | RnaOutcome
    | ProteinOutcome
    | AlleleForm
    | CodingDnaUnknown
)


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
    description: _VariantDescription


__all__ = [
    "Accession",
    "Allele",
    "AlleleStateCertainty",
    "AllelePhase",
    "AlleleVariant",
    "AlleleForm",
    "DerivedAlleleForm",
    "AlternativeAlleleForm",
    "CopiedSequence",
    "CoordinateSystem",
    "HgvsVariant",
    "Location",
    "LiteralSequence",
    "NucleotideDeletionInsertion",
    "NucleotideAnchor",
    "NucleotideCoordinate",
    "NucleotideEdit",
    "NucleotideInsertion",
    "NucleotideRepeatBlock",
    "NucleotideRepeat",
    "NucleotideSequenceItem",
    "NucleotideSubstitution",
    "NucleotideDeletion",
    "NucleotideInversion",
    "NucleotideDuplication",
    "NucleotideNoChange",
    "RnaOutcome",
    "RnaProduced",
    "RnaNoChange",
    "RnaNotProduced",
    "RnaUncertainSplicing",
    "RnaUnknown",
    "RnaIndeterminate",
    "ProteinCoordinate",
    "ProteinEdit",
    "ProteinSubstitution",
    "ProteinDeletion",
    "ProteinDeletionInsertion",
    "ProteinInsertion",
    "ProteinDuplication",
    "ProteinRepeat",
    "ProteinFrameshift",
    "ProteinExtension",
    "ProteinExtensionTerminal",
    "ProteinFrameshiftStop",
    "ProteinSequence",
    "ReferenceSpec",
    "RepeatSequenceItem",
    "VariantDescription",
    "KnownLocation",
    "KnownQuantity",
    "KnownRepeatUnit",
    "OutcomeCertainty",
    "PossibleRange",
    "Quantity",
    "Repeat",
    "RepeatUnit",
    "UncertainLocation",
    "UncertainQuantity",
    "UnknownRepeatUnit",
    "OmittedProteinFrameshiftStop",
    "UnknownProteinFrameshiftStop",
    "KnownProteinFrameshiftStop",
    "ProteinNoChange",
    "ProteinOutcome",
    "ProteinProduced",
    "ProteinNotProduced",
    "ProteinProducedAlternatives",
    "ProteinUnknown",
    "UnknownQuantity",
    "CodingDnaUnknown",
    "CodingDnaOutcome",
]
