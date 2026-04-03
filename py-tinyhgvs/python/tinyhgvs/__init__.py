"""Public Python API for the :mod:`tinyhgvs` package.

The package exports:

- :func:`parse_hgvs` as the main parsing entry point
- typed dataclasses and enums from :mod:`tinyhgvs.models`
- :class:`TinyHGVSError` and related error enums from :mod:`tinyhgvs.errors`
"""

from importlib.metadata import PackageNotFoundError, version

from .api import parse_hgvs
from .errors import ParseHgvsErrorKind, TinyHGVSError
from .models import (
    Accession,
    CopiedSequenceItem,
    CoordinateSystem,
    HgvsVariant,
    Interval,
    LiteralSequenceItem,
    NucleotideDeletionInsertionEdit,
    NucleotideAnchor,
    NucleotideCoordinate,
    NucleotideEdit,
    NucleotideInsertionEdit,
    NucleotideRepeatBlock,
    NucleotideRepeatEdit,
    NucleotideSequenceItem,
    NucleotideSequenceOmittedEdit,
    NucleotideSubstitutionEdit,
    NucleotideVariant,
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
    ReferenceSpec,
    RepeatSequenceItem,
    VariantDescription,
)

try:
    __version__ = version("tinyhgvs")
except PackageNotFoundError:
    __version__ = "0+unknown"

__all__ = [
    "Accession",
    "CopiedSequenceItem",
    "CoordinateSystem",
    "HgvsVariant",
    "Interval",
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
    "ParseHgvsErrorKind",
    "ReferenceSpec",
    "RepeatSequenceItem",
    "TinyHGVSError",
    "VariantDescription",
    "__version__",
    "parse_hgvs",
]
