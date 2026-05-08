use tinyhgvs::{
    Interval, LiteralSequenceItem, ProteinEdit, ProteinEffect, Quantity, RepeatEdit,
    RepeatSequenceUnit,
};

pub fn make_known_repeat_edit(unit: &str, count: usize) -> RepeatEdit {
    RepeatEdit {
        unit: Some(RepeatSequenceUnit::Known(LiteralSequenceItem {
            value: unit.to_string(),
        })),
        quantity: Quantity::Known { count },
    }
}

pub fn make_shorthand_repeat_edit(quantity: impl ToQuantity) -> RepeatEdit {
    RepeatEdit {
        unit: None,
        quantity: quantity.to_quantity(),
    }
}

// The "N" part inside "N[100]", "N[(100_120)]", "N[?]"
pub fn make_unknown_repeat_unit() -> RepeatSequenceUnit {
    RepeatSequenceUnit::Unknown
}

pub trait ToQuantity {
    fn to_quantity(self) -> Quantity;
}

impl ToQuantity for usize {
    fn to_quantity(self) -> Quantity {
        Quantity::Known { count: self }
    }
}

impl ToQuantity for (usize, usize) {
    fn to_quantity(self) -> Quantity {
        Quantity::Uncertain(Interval {
            start: self.0,
            end: Some(self.1),
        })
    }
}

impl ToQuantity for Quantity {
    fn to_quantity(self) -> Quantity {
        self
    }
}

pub fn make_known_quantity(count: usize) -> Quantity {
    Quantity::Known { count }
}

// The "()" part inside "N[(100_120)]" and "NM_004006.3:r.-128_-126[(600_800)]"
pub fn make_quantity_range(lo: usize, hi: usize) -> Quantity {
    Quantity::Uncertain(Interval {
        start: lo,
        end: Some(hi),
    })
}

// The "?" part inside "N[?]"
pub fn make_unknown_quantity() -> Quantity {
    Quantity::Unknown
}

pub fn get_protein_repeat(effect: &ProteinEffect) -> &RepeatEdit {
    if let ProteinEffect::Known {
        edit: ProteinEdit::Repeat(repeat),
        ..
    } = effect
    {
        repeat
    } else {
        panic!("Expected a ProteinEdit::Repeat, but found: {:?}", effect);
    }
}
