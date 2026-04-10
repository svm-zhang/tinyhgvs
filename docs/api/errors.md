# Errors

`tinyhgvs` currently exposes one public exception type, `TinyHGVSError`, with a
stable broad error kind and a more specific diagnostic code.

Representative examples:

- Invalid syntax:
  `NM_004006.2:c.5697delA` -> `invalid.syntax`
- Unsupported allele member syntax:
  `NM_004006.2:c.[2376G>C];[?]` -> `unsupported.allele_unknown_variant`
- Unsupported uncertain allele-member state:
  `NM_004006.2:c.[2376G>C](;)(1083A>C)` -> `unsupported.allele_uncertain_variant_state`
- Unsupported alternate allele-state syntax:
  `NP_003997.2:p.[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]` -> `unsupported.alternate_allele_state`
- Unsupported one-allele multi-protein syntax:
  `NP_003997.1:p.[Lys31Asn,Val25_Lys31del]` -> `unsupported.one_allele_multi_protein`
- Unsupported quantified protein insertion syntax:
  `p.Arg78_Gly79insXaa[23]` -> `unsupported.protein_insertion_payload`
- Unsupported RNA special-state syntax:
  `NM_004006.3:r.spl` -> `unsupported.rna_special_state`

::: tinyhgvs.errors
