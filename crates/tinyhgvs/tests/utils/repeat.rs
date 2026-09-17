use tinyhgvs::{
    LiteralSequenceItem, ProteinEditKind, ProteinOutcome, Quantity, RepeatEdit,
    RepeatSequenceUnit,
};

pub fn make_known_repeat_edit(unit: &str, count: usize) -> RepeatEdit {
    make_known_repeat_edit_with_quantity(unit, count)
}

pub fn make_known_repeat_edit_with_quantity(unit: &str, quantity: impl ToQuantity) -> RepeatEdit {
    RepeatEdit {
        unit: Some(RepeatSequenceUnit::Known(LiteralSequenceItem {
            value: unit.to_string(),
        })),
        quantity: quantity.to_quantity(),
    }
}

pub fn make_shorthand_repeat_edit(quantity: impl ToQuantity) -> RepeatEdit {
    RepeatEdit {
        unit: None,
        quantity: quantity.to_quantity(),
    }
}

pub fn make_unknown_repeat_edit(quantity: impl ToQuantity) -> RepeatEdit {
    RepeatEdit {
        unit: Some(RepeatSequenceUnit::Unknown),
        quantity: quantity.to_quantity(),
    }
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
        Quantity::Uncertain {
            lo: Some(self.0),
            hi: Some(self.1),
        }
    }
}

impl ToQuantity for (Option<usize>, Option<usize>) {
    fn to_quantity(self) -> Quantity {
        Quantity::Uncertain {
            lo: self.0,
            hi: self.1,
        }
    }
}

impl ToQuantity for Quantity {
    fn to_quantity(self) -> Quantity {
        self
    }
}

pub fn make_quantity_range(lo: usize, hi: usize) -> Quantity {
    Quantity::Uncertain {
        lo: Some(lo),
        hi: Some(hi),
    }
}

pub fn make_unknown_quantity() -> Quantity {
    Quantity::Unknown
}

pub fn get_protein_repeat(outcome: &ProteinOutcome) -> &RepeatEdit {
    let (edit, _) = crate::utils::ProteinOutcomeExt::produced_edit(outcome);
    if let ProteinEditKind::Repeat(repeat) = &edit.kind {
        repeat
    } else {
        panic!("expected a protein repeat");
    }
}
