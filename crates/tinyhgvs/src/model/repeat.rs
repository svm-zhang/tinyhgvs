//! Repeat units, repeat quantities, and repeat edits.

use super::LiteralSequenceItem;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepeatSequenceUnit {
    // CAG[23]
    Known(LiteralSequenceItem),
    // N[12]
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Quantity {
    // [12]
    Known {
        count: usize,
    },
    // [(60_80)], [(?_60)], [(60_?)]
    Uncertain {
        lo: Option<usize>,
        hi: Option<usize>,
    },
    // [?]
    Unknown,
}

/// Repeat edit unit and quantity.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{
///     GenomicOutcome, NucleotideEditKind, Quantity, RepeatSequenceUnit, VariantDescription,
///     parse_hgvs,
/// };
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("NC_000014.8:g.123_191CAG[19]CAA[4]")?;
///
/// let VariantDescription::Genomic(GenomicOutcome::Known(edit)) = variant.description else {
///     panic!("expected a genomic edit");
/// };
///
/// let NucleotideEditKind::Repeat { blocks } = edit.kind else {
///     panic!("expected repeat edit");
/// };
///
/// assert_eq!(blocks.len(), 2);
/// assert!(matches!(
///     &blocks[0].unit,
///     Some(RepeatSequenceUnit::Known(unit)) if unit.value == "CAG"
/// ));
/// assert_eq!(blocks[0].quantity, Quantity::Known { count: 19 });
/// # Ok(())
/// # }
/// ```
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
