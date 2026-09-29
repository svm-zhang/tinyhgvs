# Errors

`tinyhgvs` currently exposes one public exception type, `TinyHGVSError`, with a
stable broad error kind and a more specific diagnostic code.

Representative examples:

- Invalid syntax:
  `NM_004006.2:c.5697delA` -> `invalid.syntax`
- Unsupported telomeric coordinate syntax:
  `NC_000023.11:g.pter_qtersup` -> `unsupported.telomeric_position`
- Unsupported epigenetic edit syntax:
  `NC_000011.10:g.1999904_1999946|gom` -> `unsupported.epigenetic_edit`
- Unsupported RNA adjoined transcript syntax:
  `NM_002354.2:r.-358_555::NM_000251.2:r.212_*279` -> `unsupported.rna_adjoined_transcript`
- Unsupported RNA splicing outcome container:
  `NC_000023.11(NM_004006.2):r.(897u>g,832_960del)` -> `unsupported.rna_splicing_outcome`

::: tinyhgvs.errors
