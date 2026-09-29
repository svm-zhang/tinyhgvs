//! Allele models shared by genomic, coding-DNA, RNA, and protein descriptions.

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlleleStateCertainty {
    // [A] or A(;)B
    Certain,
    // A(;)(B)
    Uncertain,
}

/// Phase relationship between two established alleles.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{AlleleForm, AllelePhase, VariantDescription, parse_hgvs};
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("NM_004006.2:c.76A>G(;)80del")?;
///
/// let VariantDescription::CodingDnaAllele(AlleleForm::Single(allele)) = variant.description else {
///     panic!("expected a coding-DNA allele");
/// };
///
/// assert_eq!(allele.phase, Some(AllelePhase::Uncertain));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AllelePhase {
    // [A];[B]
    Trans,
    // A(;)B, A(;)(B)
    Uncertain,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DerivedAllele<T> {
    pub outcomes: Vec<T>,
}

impl<T> DerivedAllele<T> {
    pub fn try_from_outcomes(outcomes: Vec<T>) -> Option<Self> {
        (outcomes.len() >= 2).then_some(Self { outcomes })
    }

    pub(crate) fn from_outcomes(outcomes: Vec<T>) -> Self {
        assert!(outcomes.len() >= 2);
        Self { outcomes }
    }

    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.outcomes.iter()
    }

    pub fn map_t<U, F>(self, f: F) -> DerivedAllele<U>
    where
        F: Fn(T) -> U + Copy,
    {
        DerivedAllele {
            outcomes: self.outcomes.into_iter().map(f).collect(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlleleForm<T> {
    // [A;B], [A];[B], A(;)B
    Single(AlleleVariant<T>),
    // [A,B,C]
    Derived(DerivedAllele<T>),
    // [A]^[B]
    Alternative(Vec<AlleleVariant<T>>),
}

impl<T> AlleleForm<T> {
    pub fn map_t<U, F>(self, f: F) -> AlleleForm<U>
    where
        F: Fn(T) -> U + Copy,
    {
        match self {
            Self::Single(variant) => AlleleForm::Single(variant.map_t(f)),
            Self::Derived(derived) => AlleleForm::Derived(derived.map_t(f)),
            Self::Alternative(alternatives) => AlleleForm::Alternative(
                alternatives
                    .into_iter()
                    .map(|variant| variant.map_t(f))
                    .collect(),
            ),
        }
    }
}

/// One allele containing one or more inner variants.
///
/// Variants inside one allele are implicitly written in cis.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{AlleleForm, VariantDescription, parse_hgvs};
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("NC_000001.11:g.[123G>A;345del]")?;
///
/// let VariantDescription::GenomicAllele(AlleleForm::Single(allele)) = variant.description else {
///     panic!("expected a genomic allele");
/// };
///
/// assert_eq!(allele.allele_one.variants.len(), 2);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele<T> {
    pub variants: Vec<T>,
    pub state_certainty: AlleleStateCertainty,
}

impl<T> Allele<T> {
    /// Builds one allele from all its carrying variants.
    pub fn from_variants(variants: Vec<T>) -> Self {
        Self {
            variants,
            state_certainty: AlleleStateCertainty::Certain,
        }
    }

    pub fn uncertain_from_variants(variants: Vec<T>) -> Self {
        Self {
            variants,
            state_certainty: AlleleStateCertainty::Uncertain,
        }
    }

    /// Returns the inner variants carried by this allele.
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.variants.iter()
    }

    pub fn map_t<U, F>(self, f: F) -> Allele<U>
    where
        F: Fn(T) -> U + Copy,
    {
        Allele {
            variants: self.variants.into_iter().map(f).collect(),
            state_certainty: self.state_certainty,
        }
    }
}

impl<'a, T> IntoIterator for &'a Allele<T> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.variants.iter()
    }
}

/// Allele container holding an initial allele, an optional second established
/// allele, and any later unphased alleles.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{AlleleForm, AllelePhase, VariantDescription, parse_hgvs};
///
/// # fn main() -> Result<(), tinyhgvs::ParseHgvsError> {
/// let variant = parse_hgvs("NM_004006.2:c.[2376G>C];[2376=]")?;
///
/// let VariantDescription::CodingDnaAllele(AlleleForm::Single(allele)) = variant.description else {
///     panic!("expected a coding-DNA allele");
/// };
///
/// assert_eq!(allele.allele_one.variants.len(), 1);
/// assert!(allele.allele_two.is_some());
/// assert_eq!(allele.phase, Some(AllelePhase::Trans));
/// assert_eq!(allele.iter_outcomes().count(), 2);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleVariant<T> {
    pub allele_one: Allele<T>,
    pub allele_two: Option<Allele<T>>,
    pub phase: Option<AllelePhase>,
    pub variants_unphased: Vec<T>,
}

impl<T> AlleleVariant<T> {
    pub fn map_t<U, F>(self, f: F) -> AlleleVariant<U>
    where
        F: Fn(T) -> U + Copy,
    {
        AlleleVariant {
            allele_one: self.allele_one.map_t(f),
            allele_two: self.allele_two.map(|a| a.map_t(f)),
            phase: self.phase,
            variants_unphased: self.variants_unphased.into_iter().map(f).collect(),
        }
    }

    /// Returns all written alleles in order.
    pub fn iter_outcomes(&self) -> impl Iterator<Item = &T> {
        self.allele_one
            .variants
            .iter()
            .chain(self.allele_two.iter().flat_map(|v| v.variants.iter()))
            .chain(self.variants_unphased.iter())
    }

    /// Returns any later alleles written in uncertain relation to the
    /// established allele state.
    pub fn unphased_alleles(&self) -> &[T] {
        &self.variants_unphased
    }
}
