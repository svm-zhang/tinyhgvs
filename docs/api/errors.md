# Errors

`tinyhgvs` currently exposes one public exception type, `TinyHGVSError`, with a
stable broad error kind and a more specific diagnostic code.

Representative examples:

- Invalid syntax:
  `NM_004006.2:c.5697delA` -> `invalid.syntax`
- Unsupported allele syntax:
  `NC_000001.11:g.[123G>A;345del]` -> `unsupported.allele`
- Unsupported protein extension syntax:
  `p.(Ter157Lysext*90)` -> `unsupported.protein_extension`
- Unsupported RNA special-state syntax:
  `NM_004006.3:r.spl` -> `unsupported.rna_special_state`

::: tinyhgvs.errors
