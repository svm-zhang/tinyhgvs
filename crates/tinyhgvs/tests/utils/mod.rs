pub mod location;
pub mod repeat;

pub mod prelude {
    pub use super::location::*;
    pub use super::repeat::*;
    pub use super::{
        HgvsVariantExt, NucleotideEditKindExt, ProteinOutcomeExt, RnaOutcomeExt,
        assert_nucleotide_substitution, parse_variant,
    };
}

use tinyhgvs::{
    AlleleVariant, CodingDnaOutcome, GenomicOutcome, HgvsVariant, NucleotideEdit,
    NucleotideEditKind, NucleotideSequenceItem, OutcomeCertainty, ProteinEdit, ProteinOutcome,
    RepeatEdit, RnaOutcome, VariantDescription, parse_hgvs,
};

pub fn parse_variant(example: &str) -> HgvsVariant {
    parse_hgvs(example).unwrap_or_else(|error| panic!("{example} should parse: {error}"))
}

pub trait HgvsVariantExt {
    fn into_genomic_edit(self) -> NucleotideEdit;
    fn into_cdna_edit(self) -> NucleotideEdit;
    fn into_rna_outcome(self) -> RnaOutcome;
    fn into_protein_outcome(self) -> ProteinOutcome;
    fn into_genomic_allele(self) -> AlleleVariant<GenomicOutcome>;
    fn into_cdna_allele(self) -> AlleleVariant<CodingDnaOutcome>;
    fn into_rna_allele(self) -> AlleleVariant<RnaOutcome>;
    fn into_protein_allele(self) -> AlleleVariant<ProteinOutcome>;
}

impl HgvsVariantExt for HgvsVariant {
    fn into_genomic_edit(self) -> NucleotideEdit {
        match self.description {
            VariantDescription::Genomic(GenomicOutcome::Known(edit)) => edit,
            _ => panic!("expected a known genomic edit"),
        }
    }

    fn into_cdna_edit(self) -> NucleotideEdit {
        match self.description {
            VariantDescription::CodingDna(CodingDnaOutcome::Known(edit)) => edit,
            _ => panic!("expected a known coding-DNA edit"),
        }
    }

    fn into_rna_outcome(self) -> RnaOutcome {
        match self.description {
            VariantDescription::Rna(outcome) => outcome,
            _ => panic!("expected an RNA outcome"),
        }
    }

    fn into_protein_outcome(self) -> ProteinOutcome {
        match self.description {
            VariantDescription::Protein(outcome) => outcome,
            _ => panic!("expected a protein outcome"),
        }
    }

    fn into_genomic_allele(self) -> AlleleVariant<GenomicOutcome> {
        match self.description {
            VariantDescription::GenomicAllele(allele) => allele,
            _ => panic!("expected a genomic allele"),
        }
    }

    fn into_cdna_allele(self) -> AlleleVariant<CodingDnaOutcome> {
        match self.description {
            VariantDescription::CodingDnaAllele(allele) => allele,
            _ => panic!("expected a coding-DNA allele"),
        }
    }

    fn into_rna_allele(self) -> AlleleVariant<RnaOutcome> {
        match self.description {
            VariantDescription::RnaAllele(allele) => allele,
            _ => panic!("expected an RNA allele"),
        }
    }

    fn into_protein_allele(self) -> AlleleVariant<ProteinOutcome> {
        match self.description {
            VariantDescription::ProteinAllele(allele) => allele,
            _ => panic!("expected a protein allele"),
        }
    }
}

pub trait RnaOutcomeExt {
    fn produced_edit(&self) -> (&NucleotideEdit, &OutcomeCertainty);
}

impl RnaOutcomeExt for RnaOutcome {
    fn produced_edit(&self) -> (&NucleotideEdit, &OutcomeCertainty) {
        match self {
            Self::Produced { edit, certainty } => (edit, certainty),
            _ => panic!("expected a produced RNA outcome"),
        }
    }
}

pub trait ProteinOutcomeExt {
    fn produced_edit(&self) -> (&ProteinEdit, &OutcomeCertainty);
}

impl ProteinOutcomeExt for ProteinOutcome {
    fn produced_edit(&self) -> (&ProteinEdit, &OutcomeCertainty) {
        match self {
            Self::Produced { edit, certainty } => (edit, certainty),
            _ => panic!("expected a produced protein outcome"),
        }
    }
}

pub trait NucleotideEditKindExt {
    fn insertion_items(&self) -> &[NucleotideSequenceItem];
    fn delins_items(&self) -> &[NucleotideSequenceItem];
    fn repeat_blocks(&self) -> &[RepeatEdit];
}

impl NucleotideEditKindExt for NucleotideEditKind {
    fn insertion_items(&self) -> &[NucleotideSequenceItem] {
        match self {
            Self::Insertion { items } => items,
            _ => panic!("expected a nucleotide insertion"),
        }
    }

    fn delins_items(&self) -> &[NucleotideSequenceItem] {
        match self {
            Self::DeletionInsertion { items } => items,
            _ => panic!("expected a nucleotide deletion-insertion"),
        }
    }

    fn repeat_blocks(&self) -> &[RepeatEdit] {
        match self {
            Self::Repeat { blocks } => blocks,
            _ => panic!("expected a nucleotide repeat"),
        }
    }
}

pub fn assert_nucleotide_substitution(edit: &NucleotideEdit, reference: &str, alternate: &str) {
    assert!(matches!(
        &edit.kind,
        NucleotideEditKind::Substitution {
            reference: observed_reference,
            alternate: observed_alternate,
        } if observed_reference == reference && observed_alternate == alternate
    ));
}
