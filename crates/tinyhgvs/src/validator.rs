use crate::model::{
    CoordinateSystem, Location, NucleotideCoordinate, NucleotideEdit, NucleotideVariant,
    RepeatEdit, RepeatSequenceUnit, VariantDescription,
};

// Will likely refactor these 2 helpers to the validation module
fn is_invalid_genomic_repeat(block: &RepeatEdit) -> bool {
    !block.is_unit_known() || !block.is_copy_known()
}

fn is_invalid_cdna_repeat(block: &RepeatEdit) -> bool {
    if is_invalid_genomic_repeat(block) {
        return true;
    }

    // Invalid: NM_024312.4:c.2686A[10]
    // Invalid: NM_024312.4:c.1738TA[6]
    if let Some(RepeatSequenceUnit::Known(ref item)) = block.unit {
        return item.value.len() % 3 != 0;
    }
    false
}

fn is_invalid_rna_repeat(
    location: &Location<NucleotideCoordinate>,
    repeat_blocks: &[RepeatEdit],
) -> bool {
    // Invalid: r.-125_-123cug[4]
    // Valid: r.-125_-123[4]
    if repeat_blocks.len() == 1 {
        let known_interval_location = location.is_interval() && !location.is_uncertain();
        let is_repeat_unit_known = repeat_blocks[0].is_unit_known();
        if known_interval_location && is_repeat_unit_known {
            return true;
        }
        return false;
    }
    // https://github.com/HGVSnomenclature/hgvs-nomenclature/issues/114
    // Invalid: r.456_465[4]466_489[9]490_499[3]
    // Valid: r.456_499us[4]cag[9]gccag[3]
    repeat_blocks.iter().any(|block| !block.is_unit_known())
}

fn validate_nucleotide_repeat(
    coordinate_system: CoordinateSystem,
    location: &Location<NucleotideCoordinate>,
    repeat_blocks: &[RepeatEdit],
) -> bool {
    match coordinate_system {
        CoordinateSystem::Rna => is_invalid_rna_repeat(location, repeat_blocks),
        // Valid: NM_002024.5:c.-128_-69GGC[10]GGA[1]GGC[9]GGA[1]GGC[10]
        CoordinateSystem::CodingDna => repeat_blocks.iter().any(is_invalid_cdna_repeat),
        // Valid: NC_000012.11:g.112036755_112036823CTG[9]TTG[1]CTG[13]
        // Valid: NC_000014.8:g.101179660_101179695TG[14]
        CoordinateSystem::Genomic => repeat_blocks.iter().any(is_invalid_genomic_repeat),
        _ => false,
    }
}

fn validate_nucleotide_variant_description(
    coordinate_system: CoordinateSystem,
    location: &Location<NucleotideCoordinate>,
    edit: &NucleotideEdit,
) -> bool {
    match edit {
        NucleotideEdit::Repeat { blocks } => {
            validate_nucleotide_repeat(coordinate_system, location, blocks)
        }
        _ => false,
    }
}

fn validate_nucleotide_allele_description(
    coordinate_system: CoordinateSystem,
    variants: &[NucleotideVariant],
) -> bool {
    match coordinate_system {
        CoordinateSystem::Genomic => variants.iter().any(|v| match &v.edit {
            NucleotideEdit::Repeat { blocks } => blocks.iter().any(is_invalid_genomic_repeat),
            _ => false,
        }),
        _ => false,
    }
}

pub fn validate_nucleotide_description(
    coordinate_system: CoordinateSystem,
    description: &VariantDescription,
) -> bool {
    match description {
        VariantDescription::Nucleotide(NucleotideVariant { location, edit }) => {
            validate_nucleotide_variant_description(coordinate_system, location, edit)
        }
        VariantDescription::NucleotideAllele(allele) => {
            validate_nucleotide_allele_description(coordinate_system, &allele.allele_one.variants)
                || allele.allele_two.as_ref().map_or(false, |allele| {
                    validate_nucleotide_allele_description(coordinate_system, &allele.variants)
                })
        }
        _ => false,
    }
}
