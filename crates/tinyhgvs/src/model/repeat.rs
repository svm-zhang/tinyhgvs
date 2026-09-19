use super::LiteralSequenceItem;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepeatSequenceUnit {
    Known(LiteralSequenceItem), // reuse the LiteralSequenceItem struct
    Unknown,                    // N or n
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Quantity {
    Known {
        count: usize,
    },
    Uncertain {
        lo: Option<usize>,
        hi: Option<usize>,
    },
    // Uncertain(Interval<usize>), // reuse the Interval type with usize.
    Unknown, // [?] case
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatEdit {
    pub quantity: Quantity,
    pub unit: Option<RepeatSequenceUnit>,
}

impl RepeatEdit {
    pub fn is_unit_known(&self) -> bool {
        matches!(self.unit, Some(RepeatSequenceUnit::Known(_)))
    }

    pub fn is_copy_known(&self) -> bool {
        matches!(self.quantity, Quantity::Known { .. })
    }

    pub fn is_copy_unknown(&self) -> bool {
        matches!(self.quantity, Quantity::Unknown)
    }
}
