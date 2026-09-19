use super::NucleotideCoordinate;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval<T> {
    pub start: T,
    pub end: Option<T>,
}

impl Interval<NucleotideCoordinate> {
    fn is_end_bound_unknown(&self) -> bool {
        self.end
            .as_ref()
            .map_or(false, NucleotideCoordinate::is_unknown)
    }

    /// Returns `true` when either bound of the interval is unknown, i.e. `?_B`,
    /// `A_?`, `?_?`
    pub fn has_unknown_bound(&self) -> bool {
        self.start.is_unknown() || self.is_end_bound_unknown()
    }

    /// Returns `true` when both sides of the interval are unknown, i.e. ?_?
    pub fn is_fully_unknown(&self) -> bool {
        self.start.is_unknown() && self.is_end_bound_unknown()
    }
}

/// Main edited location on a nucleotide or protein variant/effect.
///
/// Known locations keep the current one-level interval shape. Uncertain
/// locations wrap the left and right uncertain regions as intervals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Location<T> {
    Known(Interval<T>),
    Uncertain(Interval<Interval<T>>),
}

impl<T> Location<T> {
    /// Builds one known location from a plain interval.
    pub fn from_known(value: Interval<T>) -> Self {
        Self::Known(value)
    }

    /// Builds one uncertain location from uncertain left/right regions.
    pub fn from_uncertain(value: Interval<Interval<T>>) -> Self {
        Self::Uncertain(value)
    }

    /// Returns `true` when the location is written with uncertain-region
    /// syntax.
    pub fn is_uncertain(&self) -> bool {
        matches!(self, Self::Uncertain(_))
    }

    /// Returns `true` for one known position.
    pub fn is_pos(&self) -> bool {
        matches!(self, Self::Known(interval) if interval.end.is_none())
    }

    /// Returns `true` for interval-shaped locations, whether known or
    /// uncertain.
    pub fn is_interval(&self) -> bool {
        !self.is_pos()
    }

    /// Returns the left known position when the location is known.
    pub fn start(&self) -> Option<&T> {
        match self {
            Self::Known(interval) => Some(&interval.start),
            Self::Uncertain(_) => None,
        }
    }

    /// Returns the right known position when the location is a known interval.
    pub fn end(&self) -> Option<&T> {
        match self {
            Self::Known(interval) => interval.end.as_ref(),
            Self::Uncertain(_) => None,
        }
    }

    /// Returns the left uncertain region when the location is uncertain.
    pub fn l_interval(&self) -> Option<&Interval<T>> {
        match self {
            Self::Known(_) => None,
            Self::Uncertain(interval) => Some(&interval.start),
        }
    }

    /// Returns the right uncertain region when the location is an uncertain
    /// interval.
    pub fn r_interval(&self) -> Option<&Interval<T>> {
        match self {
            Self::Known(_) => None,
            Self::Uncertain(interval) => interval.end.as_ref(),
        }
    }
}
