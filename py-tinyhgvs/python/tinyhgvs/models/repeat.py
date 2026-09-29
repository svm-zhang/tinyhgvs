from dataclasses import dataclass
from typing import TypeAlias


@dataclass(frozen=True, slots=True)
class KnownQuantity:
    """Known repeat quantity.

    Examples:
        >>> from tinyhgvs import KnownQuantity, parse_hgvs
        >>> variant = parse_hgvs("NC_000014.8:g.123CAG[23]")
        >>> quantity = variant.description.sequence[0].quantity
        >>> isinstance(quantity, KnownQuantity)
        True
        >>> quantity.count
        23
    """

    count: int


@dataclass(frozen=True, slots=True)
class UncertainQuantity:
    """Uncertain repeat quantity range.

    Examples:
        >>> from tinyhgvs import UncertainQuantity, parse_hgvs
        >>> variant = parse_hgvs("NC_000003.12:g.63912687AGC[(19_23)]")
        >>> quantity = variant.description.sequence[0].quantity
        >>> isinstance(quantity, UncertainQuantity)
        True
        >>> quantity.lo, quantity.hi
        (19, 23)
    """

    lo: int | None
    hi: int | None


@dataclass(frozen=True, slots=True)
class UnknownQuantity:
    """Unknown repeat quantity written as ``[?]``.

    Examples:
        >>> from tinyhgvs import UnknownQuantity, parse_hgvs
        >>> variant = parse_hgvs("NC_000003.12:g.63912687AGC[?]")
        >>> quantity = variant.description.sequence[0].quantity
        >>> isinstance(quantity, UnknownQuantity)
        True
    """

    pass


Quantity: TypeAlias = KnownQuantity | UncertainQuantity | UnknownQuantity
"""Tagged union for supported repeat quantity models:

- [`KnownQuantity`][tinyhgvs.models.repeat.KnownQuantity]
- [`UncertainQuantity`][tinyhgvs.models.repeat.UncertainQuantity]
- [`UnknownQuantity`][tinyhgvs.models.repeat.UnknownQuantity]
"""


@dataclass(frozen=True, slots=True)
class KnownRepeatUnit:
    """Known repeated sequence unit.

    Examples:
        >>> from tinyhgvs import KnownRepeatUnit, parse_hgvs
        >>> variant = parse_hgvs("NC_000014.8:g.123CAG[23]")
        >>> unit = variant.description.sequence[0].unit
        >>> isinstance(unit, KnownRepeatUnit)
        True
        >>> unit.value
        'CAG'
    """

    value: str


@dataclass(frozen=True, slots=True)
class UnknownRepeatUnit:
    """Unknown repeated sequence unit.

    Examples:
        >>> from tinyhgvs import UnknownRepeatUnit, parse_hgvs
        >>> variant = parse_hgvs("NC_000023.10:g.32717298_32717299insN[100]")
        >>> unit = variant.description.sequence[0].unit
        >>> isinstance(unit, UnknownRepeatUnit)
        True
    """

    pass


RepeatUnit: TypeAlias = KnownRepeatUnit | UnknownRepeatUnit
"""Tagged union for supported repeat unit models:

- [`KnownRepeatUnit`][tinyhgvs.models.repeat.KnownRepeatUnit]
- [`UnknownRepeatUnit`][tinyhgvs.models.repeat.UnknownRepeatUnit]
"""


@dataclass(frozen=True, slots=True)
class Repeat:
    """Repeat unit and quantity used by nucleotide and protein repeat syntax.

    A repeat unit can be known, unknown, or omitted when the unit is implied by
    the surrounding location.

    Examples:
        A literal nucleotide repeat unit:
        >>> from tinyhgvs import parse_hgvs
        >>> variant = parse_hgvs("NC_000014.8:g.123CAG[23]")
        >>> repeat = variant.description.sequence[0]
        >>> repeat.unit.value
        'CAG'
        >>> repeat.quantity.count
        23
        >>> repeat.is_unit_known
        True

        A repeat count can be unknown:
        >>> variant = parse_hgvs("NC_000003.12:g.63912687AGC[?]")
        >>> repeat = variant.description.sequence[0]
        >>> repeat.is_quantity_unknown
        True

        A repeat count can be uncertain:
        >>> variant = parse_hgvs("NC_000003.12:g.63912687AGC[(19_23)]")
        >>> repeat = variant.description.sequence[0]
        >>> repeat.is_quantity_uncertain
        True

        Protein repeats use the same repeat model:
        >>> variant = parse_hgvs("NP_0123456.1:p.Arg65_Ser67[12]")
        >>> variant.description.edit.repeat.quantity.count
        12
    """

    unit: RepeatUnit | None
    quantity: Quantity

    @property
    def is_unit_known(self) -> bool:
        return isinstance(self.unit, KnownRepeatUnit)

    @property
    def is_quantity_known(self) -> bool:
        return isinstance(self.quantity, KnownQuantity)

    @property
    def is_quantity_uncertain(self) -> bool:
        return isinstance(self.quantity, UncertainQuantity)

    @property
    def is_quantity_unknown(self) -> bool:
        return isinstance(self.quantity, UnknownQuantity)


__all__ = [
    "KnownQuantity",
    "Quantity",
    "UncertainQuantity",
    "UnknownQuantity",
    "Repeat",
    "KnownRepeatUnit",
    "UnknownRepeatUnit",
    "RepeatUnit",
]
