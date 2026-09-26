from dataclasses import dataclass

from .nucleotide import NucleotideEdit
from .shared import OutcomeCertainty


class RnaOutcome:
    """Base class for RNA outcomes."""

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
    def is_uncertain_splicing(self) -> bool:
        return False

    @property
    def is_unknown(self) -> bool:
        return False

    @property
    def is_indeterminate(self) -> bool:
        return False


@dataclass(frozen=True, slots=True)
class RnaProduced(RnaOutcome):
    """RNA outcome where a nucleotide edit is produced.

    Examples:
        >>> from tinyhgvs import RnaProduced, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.3:r.456_465del")
        >>> isinstance(variant.description, RnaProduced)
        True
        >>> variant.description.edit.location.start.coordinate
        456
    """

    edit: NucleotideEdit
    certainty: OutcomeCertainty

    @property
    def is_produced(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED


@dataclass(frozen=True, slots=True)
class RnaNoChange(RnaOutcome):
    """RNA outcome with no sequence change.

    Examples:
        >>> from tinyhgvs import RnaNoChange, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.3:r.=")
        >>> isinstance(variant.description, RnaNoChange)
        True
        >>> variant.description.is_no_change
        True
    """

    certainty: OutcomeCertainty

    @property
    def is_no_change(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED


@dataclass(frozen=True, slots=True)
class RnaNotProduced(RnaOutcome):
    """RNA outcome where no RNA product is produced.

    Examples:
        >>> from tinyhgvs import RnaNotProduced, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.3:r.0")
        >>> isinstance(variant.description, RnaNotProduced)
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
class RnaUncertainSplicing(RnaOutcome):
    """RNA outcome for uncertain splicing, written as ``r.spl``.

    Examples:
        >>> from tinyhgvs import RnaUncertainSplicing, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.3:r.spl")
        >>> isinstance(variant.description, RnaUncertainSplicing)
        True
    """

    @property
    def is_uncertain_splicing(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class RnaUnknown(RnaOutcome):
    """RNA outcome ``r.?`` where the RNA consequence is unknown.

    Examples:
        >>> from tinyhgvs import RnaUnknown, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.3:r.?")
        >>> isinstance(variant.description, RnaUnknown)
        True
    """

    @property
    def is_unknown(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class RnaIndeterminate(RnaOutcome):
    """RNA outcome ``r.(?)`` where the RNA consequence is indeterminate.

    Examples:
        >>> from tinyhgvs import RnaIndeterminate, parse_hgvs
        >>> variant = parse_hgvs("NM_004006.3:r.(?)")
        >>> isinstance(variant.description, RnaIndeterminate)
        True
    """

    @property
    def is_indeterminate(self) -> bool:
        return True


__all__ = [
    "RnaOutcome",
    "RnaProduced",
    "RnaNoChange",
    "RnaNotProduced",
    "RnaUncertainSplicing",
    "RnaUnknown",
    "RnaIndeterminate",
]
