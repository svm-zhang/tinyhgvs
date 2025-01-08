use pest::Parser;
use pest_derive::Parser;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[derive(Parser)]
#[grammar = "grammar.pest"]
pub struct HgvsParser;

#[pyclass]
#[derive(Debug, Clone)]
pub struct HgvsPosition {
    #[pyo3(get)]
    pub prefix: String,
    #[pyo3(get)]
    pub base: i32,
    #[pyo3(get)]
    pub offset: i32,
}

impl HgvsPosition {
    pub fn new(prefix: Option<String>, base: i32, offset: Option<i32>) -> Self {
        HgvsPosition {
            prefix: prefix.unwrap_or_default(),
            base,
            offset: offset.unwrap_or(0),
        }
    }
}

#[pyclass]
#[derive(Debug, PartialEq, Clone)]
pub struct HgvsEdit {
    #[pyo3(get)]
    pub edit_type: String,
    #[pyo3(get)]
    pub reference: String,
    #[pyo3(get)]
    pub alt: String,
}

impl HgvsEdit {
    pub fn new() -> Self {
        HgvsEdit {
            edit_type: String::new(),
            reference: String::new(),
            alt: String::new(),
        }
    }
}

impl Default for HgvsEdit {
    fn default() -> Self {
        Self::new()
    }
}

#[pyclass]
#[derive(Debug, Clone)]
pub struct HgvsVariant {
    #[pyo3(get)]
    pub ac: String,
    #[pyo3(get)]
    pub seq_type: String,
    #[pyo3(get)]
    pub start: HgvsPosition,
    #[pyo3(get)]
    pub end: Option<HgvsPosition>,
    #[pyo3(get)]
    pub edit: HgvsEdit,
}

fn parse_hgvsd() -> bool {
    !todo!()
}

#[pyfunction]
pub fn parse_hgvs(input: &str) -> PyResult<HgvsVariant> {
    let pairs = HgvsParser::parse(Rule::hgvsc, input)
        .map_err(|_| PyValueError::new_err("Failed to parse rule"))?
        .next()
        .unwrap();
    let mut ac = String::new();
    let mut seq_type = String::new();
    let mut pos = HgvsPosition::new(None, 0, None);
    let mut end = HgvsPosition::new(None, 0, None);
    let mut edit = HgvsEdit::default();

    for pair in pairs.into_inner() {
        match pair.as_rule() {
            Rule::ac => ac = pair.as_str().to_string(),
            Rule::coord_type => seq_type = pair.as_str().to_string(),
            Rule::pos => {
                let mut prefix = None;
                let mut base = 0;
                let mut offset = None;
                let mut offset_direction = 0;
                for sub_pair in pair.into_inner() {
                    match sub_pair.as_rule() {
                        Rule::utr_prefix => {
                            prefix = Some(sub_pair.as_str().to_string());
                            offset_direction = if sub_pair.as_str() == "-" { -1 } else { 1 };
                        }
                        Rule::base => {
                            base = sub_pair.as_str().parse::<i32>().unwrap();
                        }
                        Rule::offset => {
                            offset = Some(sub_pair.as_str().parse::<i32>().unwrap());
                            if offset_direction != 0 {
                                offset = Some(offset.unwrap() * offset_direction);
                            }
                        }
                        Rule::intronic_prefix => {
                            offset_direction = if sub_pair.as_str() == "-" { -1 } else { 1 };
                        }
                        _ => {}
                    }
                }
                pos = HgvsPosition::new(prefix.clone(), base, offset);
                // end.prefix = prefix.clone().unwrap();
                // end.base = base;
                // end.offset = offset.unwrap();
            }
            Rule::range => {
                let mut start_base = 0;
                let mut end_base = 0;
                let mut offset = None;
                let mut offset_direction = 0;
                for sub_pair in pair.into_inner() {
                    match sub_pair.as_rule() {
                        Rule::range_nuc => {
                            let mut found_parts = sub_pair.into_inner();
                            start_base = found_parts
                                .next()
                                .map_or(0, |part| part.as_str().parse::<i32>().unwrap());
                            println!("{}", start_base);
                            end_base = found_parts
                                .next()
                                .map_or(0, |part| part.as_str().parse::<i32>().unwrap());
                        }
                        Rule::intronic_prefix => {
                            offset_direction = if sub_pair.as_str() == "-" { -1 } else { 1 };
                        }
                        Rule::offset => {
                            offset = Some(sub_pair.as_str().parse::<i32>().unwrap());
                            if offset_direction != 0 {
                                offset = Some(offset.unwrap() * offset_direction);
                            }
                        }
                        _ => {}
                    }
                }
                // pos = HgvsPosition::new(None, start_base, Some(0));
                // end = HgvsPosition::new(None, end_base, offset)
                pos.base = start_base;
                end.base = end_base;
                end.offset = offset.unwrap_or(0);
            }
            Rule::sub => {
                let mut found_parts = pair.into_inner();
                edit.edit_type = "sub".to_string();
                edit.reference = found_parts.next().unwrap().as_str().to_string();
                edit.alt = found_parts.next().unwrap().as_str().to_string();
            }
            Rule::deletion => {
                let mut found_parts = pair.into_inner();
                edit.edit_type = found_parts.next().unwrap().as_str().to_string();
                edit.reference = String::new();
                edit.alt = found_parts
                    .next()
                    .map_or(String::new(), |part| part.as_str().to_string());
            }
            _ => {}
        }
    }
    Ok(HgvsVariant {
        ac,
        seq_type,
        start: pos,
        end: Some(end),
        edit,
    })
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_exonic_substitution() {
        let hgvs = "NM_000001:c.1234A>G";
        let variant = parse_hgvs(hgvs).unwrap();
        assert_eq!(variant.ac, "NM_000001");
        assert_eq!(variant.seq_type, "c");
        assert_eq!(variant.start.base, 1234);
        assert_eq!(variant.start.offset, 0);
        assert_eq!(variant.edit.reference, "A");
        assert_eq!(variant.edit.alt, "G");
    }

    #[test]
    fn test_following_intronic_substitution() {
        let hgvs = "NM_004006:c.93+1C>T";
        let variant = parse_hgvs(hgvs).unwrap();
        println!(
            "{} {} {} {}",
            hgvs, variant.start.prefix, variant.start.base, variant.start.offset
        );
        assert_eq!(variant.ac, "NM_004006");
        assert_eq!(variant.seq_type, "c");
        assert_eq!(variant.start.base, 93);
        assert_eq!(variant.start.offset, 1);
        assert_eq!(variant.edit.reference, "C");
        assert_eq!(variant.edit.alt, "T");
    }

    #[test]
    fn test_preceding_intronic_substitution() {
        let hgvs = "NM_000002:c.371-18C>T";
        let variant = parse_hgvs(hgvs).unwrap();
        assert_eq!(variant.ac, "NM_000002");
        assert_eq!(variant.seq_type, "c");
        assert_eq!(variant.start.base, 371);
        assert_eq!(variant.start.offset, -18);
    }

    #[test]
    fn test_utr5_substitution() {
        // substitution in the 5' prime UTR
        let hgvs = "NM_000003:c.-18A>T";
        let variant = parse_hgvs(hgvs).unwrap();
        println!(
            "{} {} {}",
            variant.start.prefix, variant.start.base, variant.start.offset
        );
        assert_eq!(variant.ac, "NM_000003");
        assert_eq!(variant.seq_type, "c");
        assert_eq!(variant.start.prefix, "-");
        assert_eq!(variant.start.base, 0);
        assert_eq!(variant.start.offset, -18);
    }

    #[test]
    fn test_utr3_substitution() {
        // substitution in the 3' prime UTR
        let hgvs = "NM_000003:c.*18A>T";
        let variant = parse_hgvs(hgvs).unwrap();
        assert_eq!(variant.ac, "NM_000003");
        assert_eq!(variant.seq_type, "c");
        assert_eq!(variant.start.prefix, "*");
        assert_eq!(variant.start.base, 0);
        assert_eq!(variant.start.offset, 18);
    }

    #[test]
    fn test_basic_one_bp_deletion_without_alt_seq() {
        let hgvs = "NM_004006:c.5697del";
        let variant = parse_hgvs(hgvs).unwrap();
        println!(
            "{} {} {} {}",
            hgvs, variant.edit.edit_type, variant.edit.reference, variant.edit.alt,
        );
        assert_eq!(variant.start.prefix, "");
        assert_eq!(variant.start.base, 5697);
        assert_eq!(variant.start.offset, 0);
        assert_eq!(variant.edit.edit_type, "del");
        assert_eq!(variant.edit.reference, "");
        assert_eq!(variant.edit.alt, "");
    }

    #[test]
    fn test_basic_one_bp_deletion_with_alt_seq() {
        let hgvs = "NM_004006:c.5697delA";
        let variant = parse_hgvs(hgvs).unwrap();
        println!(
            "{} {} {} {}",
            hgvs, variant.edit.edit_type, variant.edit.reference, variant.edit.alt,
        );
        assert_eq!(variant.start.prefix, "");
        assert_eq!(variant.start.base, 5697);
        assert_eq!(variant.start.offset, 0);
        assert_eq!(variant.edit.edit_type, "del");
        assert_eq!(variant.edit.reference, "");
        assert_eq!(variant.edit.alt, "A");
    }

    #[test]
    fn test_multi_bp_exonic_deletion_without_alt_seq() {
        let hgvs = "NM_015308:c.172_177del";
        let variant = parse_hgvs(hgvs).unwrap();
        // if let Ok(variant) = parse_hgvs(hgvs) {
        //     println!("success");
        //     println!("{} {} {}", hgvs, variant.start.base, variant.start.offset);
        //     println!("{:?}", variant.end);
        //     if variant.end.is_some() {
        //         let end = variant.end.clone().unwrap();
        //         println!("{} {}", end.base, end.offset)
        //     }
        // } else {
        //     println!("Failed parsing {}", hgvs);
        // }
        assert_eq!(variant.start.base, 172);
        assert!(variant.end.is_some());
        assert_eq!(variant.end.unwrap().base, 177);
        assert_eq!(variant.edit.edit_type, "del");
    }

    #[test]
    fn test_multi_bp_intronic_deletion_without_alt_seq() {
        let hgvs = "NM_004006:c.183_186+48del";
        let variant = parse_hgvs(hgvs).unwrap();
        if let Ok(variant) = parse_hgvs(hgvs) {
            println!("success");
            println!("{} {} {}", hgvs, variant.start.base, variant.start.offset);
            println!("{:?}", variant.end);
            if variant.end.is_some() {
                let end = variant.end.clone().unwrap();
                println!("{} {}", end.base, end.offset)
            }
        } else {
            println!("Failed parsing {}", hgvs);
        }
        let end = variant.end.unwrap();
        assert_eq!(variant.start.base, 183);
        assert_eq!(end.base, 186);
        assert_eq!(variant.start.offset, 0);
        assert_eq!(end.offset, 48);
        assert_eq!(variant.edit.edit_type, "del");
    }

    //
    // #[test]
    // fn fuck() {
    //     let input = "NM_015308:c.172_177del";
    //     let pairs = HgvsParser::parse(Rule::hgvs, input).expect("Failed to parse");
    //     for pair in pairs {
    //         println!("Rule: {:?}", pair.as_rule());
    //         println!("texg: {:?}", pair.as_str());
    //     }
    // }

    #[test]
    fn test_one_bp_intronic_deletion_without_alt_seq() {
        let hgvs = "LRG_199t1:c.1704+1del";
        let variant = parse_hgvs(hgvs).unwrap();
        println!(
            "{} {} {} {}",
            hgvs, variant.edit.edit_type, variant.edit.reference, variant.edit.alt,
        );
        assert_eq!(variant.ac, "LRG_199t1");
        assert_eq!(variant.start.prefix, "");
        assert_eq!(variant.start.base, 1704);
        assert_eq!(variant.start.offset, 1);
        assert_eq!(variant.edit.edit_type, "del");
        assert_eq!(variant.edit.reference, "");
        assert_eq!(variant.edit.alt, "");
    }
}
