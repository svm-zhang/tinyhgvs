from dataclasses import dataclass
from typing import TypeAlias


@dataclass(frozen=True, slots=True)
class KnownQuantity:
    count: int


@dataclass(frozen=True, slots=True)
class UncertainQuantity:
    lo: int | None
    hi: int | None


@dataclass(frozen=True, slots=True)
class UnknownQuantity:
    pass


Quantity: TypeAlias = KnownQuantity | UncertainQuantity | UnknownQuantity


@dataclass(frozen=True, slots=True)
class KnownRepeatUnit:
    value: str


@dataclass(frozen=True, slots=True)
class UnknownRepeatUnit:
    pass


RepeatUnit: TypeAlias = KnownRepeatUnit | UnknownRepeatUnit


@dataclass(frozen=True, slots=True)
class Repeat:
    unit: RepeatUnit | None
    quantity: Quantity

    @property
    def is_unit_known(self) -> bool:
        return isinstance(self.unit, KnownRepeatUnit)

    @property
    def is_copy_known(self) -> bool:
        return isinstance(self.quantity, KnownQuantity)

    @property
    def is_copy_unknown(self) -> bool:
        return isinstance(self.quantity, UnknownQuantity)


__all__ = [
    "KnownQuantity",
    "UncertainQuantity",
    "UnknownQuantity",
    "Repeat",
    "KnownRepeatUnit",
    "UnknownRepeatUnit",
    "RepeatUnit",
]
