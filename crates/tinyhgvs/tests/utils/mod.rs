pub mod repeat;

pub mod prelude {
    // Re-export the traits so the methods (.into_protein_variant()) work
    pub use super::VariantDescriptionExt;

    // Re-export all the helper functions from repeat
    pub use super::repeat::*;
}

use tinyhgvs::{HgvsVariant, NucleotideVariant, ProteinVariant, VariantDescription};

pub trait VariantDescriptionExt {
    fn into_nucleotide_variant(self) -> NucleotideVariant;
    fn into_protein_variant(self) -> ProteinVariant;
}

impl VariantDescriptionExt for HgvsVariant {
    fn into_protein_variant(self) -> ProteinVariant {
        match self.description {
            VariantDescription::Protein(p) => p,
            _ => panic!("Not a protein variant"),
        }
    }

    fn into_nucleotide_variant(self) -> NucleotideVariant {
        match self.description {
            VariantDescription::Nucleotide(n) => n,
            _ => panic!("Not a nucleotide variant"),
        }
    }
}
