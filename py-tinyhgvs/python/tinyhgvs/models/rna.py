from dataclasses import dataclass

from .nucleotide import NucleotideEdit
from .shared import OutcomeCertainty


class RnaOutcome:
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
    certainty: OutcomeCertainty

    @property
    def is_no_change(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED


@dataclass(frozen=True, slots=True)
class RnaNotProduced(RnaOutcome):
    certainty: OutcomeCertainty

    @property
    def is_not_produced(self) -> bool:
        return True

    @property
    def is_predicted(self) -> bool:
        return self.certainty is OutcomeCertainty.PREDICTED


@dataclass(frozen=True, slots=True)
class RnaUncertainSplicing(RnaOutcome):
    @property
    def is_uncertain_splicing(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class RnaUnknown(RnaOutcome):
    @property
    def is_unknown(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class RnaIndeterminate(RnaOutcome):
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
