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
- Unsupported protein allele syntax:
  `NP_003997.1:p.Val7=/del` -> `unsupported.protein_allele`
- Unsupported quantified protein insertion syntax:
  `p.Arg78_Gly79insXaa[23]` -> `unsupported.protein_insertion_payload`
- Unsupported RNA special-state syntax:
  `NM_004006.3:r.spl` -> `unsupported.rna_special_state`

::: tinyhgvs.errors
