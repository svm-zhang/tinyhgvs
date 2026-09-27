# Unsupported Syntax

`tinyhgvs` parses supported HGVS syntax into structured Rust and Python models.
For a small set of recognized HGVS families that are not modeled yet, parsing
raises a structured unsupported-syntax error instead of a generic syntax error.

Unsupported diagnostic codes begin with `unsupported.`. The current active
unsupported families are listed below.

## What Error To Expect

Python:

```python
from tinyhgvs import TinyHGVSError, parse_hgvs

try:
    parse_hgvs("NC_000023.11:g.pter_qtersup")
except TinyHGVSError as error:
    print(error.kind.value)
    print(error.code)
    print(error.fragment)
```

Expected output:

```text
unsupported_syntax
unsupported.telomeric_position
pter
```

Rust:

```rust
let error = tinyhgvs::parse_hgvs("NC_000023.11:g.pter_qtersup").unwrap_err();
assert_eq!(error.code(), "unsupported.telomeric_position");
assert_eq!(error.fragment(), Some("pter"));
```

## `unsupported.telomeric_position`

Telomeric DNA coordinate syntax is recognized but not modeled yet.

Examples:

- `NC_000023.11:g.pter_qtersup`
- `NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]`

Expected error:

```text
kind: unsupported_syntax
code: unsupported.telomeric_position
fragment: pter or qter
```

## `unsupported.epigenetic_edit`

DNA epigenetic edit modifiers are recognized but not modeled yet.

Example:

- `NC_000011.10:g.1999904_1999946|gom`

Expected error:

```text
kind: unsupported_syntax
code: unsupported.epigenetic_edit
fragment: |gom
```

## `unsupported.rna_adjoined_transcript`

RNA adjoined transcript syntax is recognized but not modeled yet.

Examples:

- `NM_002354.2:r.-358_555::NM_000251.2:r.212_*279`
- `NM_152263.2:r.-115_775::aggcucccuugg::NM_002609.3:r.1580_*1924`

Expected error:

```text
kind: unsupported_syntax
code: unsupported.rna_adjoined_transcript
fragment: ::
```

## `unsupported.rna_splicing_outcome`

Some higher-level RNA splicing outcome containers remain unsupported. Direct
RNA outcomes such as `r.?`, `r.(?)`, `r.0`, `r.0?`, and `r.spl` are supported.
Derived RNA allele forms such as `r.[897u>g,832_960del]` are also supported.

Example:

- `NC_000023.11(NM_004006.2):r.(897u>g,832_960del)`

Expected error:

```text
kind: unsupported_syntax
code: unsupported.rna_splicing_outcome
fragment: r.(...)
```

