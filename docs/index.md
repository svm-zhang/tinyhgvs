# tinyhgvs

[![PyPI](https://img.shields.io/pypi/v/tinyhgvs.svg)](https://pypi.org/project/tinyhgvs/)
[![Python Versions](https://img.shields.io/pypi/pyversions/tinyhgvs.svg)](https://pypi.org/project/tinyhgvs/)
[![CI](https://github.com/svm-zhang/tinyhgvs/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/svm-zhang/tinyhgvs/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/svm-zhang/tinyhgvs.svg)](https://github.com/svm-zhang/tinyhgvs/blob/main/LICENSE)
[![crates.io](https://img.shields.io/crates/v/tinyhgvs.svg)](https://crates.io/crates/tinyhgvs)
[![docs.rs](https://img.shields.io/docsrs/tinyhgvs)](https://docs.rs/tinyhgvs)

`tinyhgvs` is a lightweight HGVS parsing library with a Rust core and a Python
API. The project focuses on a small but expandable data model, clear error
diagnostics for unsupported syntax, and parity between the Rust and Python
surfaces.

## Installation

`tinyhgvs` can be installed directly from PyPI:

```bash
pip install tinyhgvs

```

### Local deployment

For local development of the Python package:

```bash
cd py-tinyhgvs
uv run maturin develop
```


## Quick Examples

Parse a coding DNA substitution:

```python
from tinyhgvs import parse_hgvs

variant = parse_hgvs("NM_004006.2:c.357+1G>A")
print(variant.coordinate_system.value)
print(variant.description.edit)
```

Parse an exact repeat:

```python
from tinyhgvs import parse_hgvs

variant = parse_hgvs("NM_004006.3:r.-124_-123[14]")
print(variant.description.edit.blocks[0].count)
```

Parse a protein frameshift consequence:

```python
from tinyhgvs import parse_hgvs

variant = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
print(variant.description.effect.edit.to_residue)
print(variant.description.effect.edit.stop.ordinal)
```

Inspect an unsupported syntax error:

```python
from tinyhgvs import TinyHGVSError, parse_hgvs

try:
    parse_hgvs("NP_003997.1:p.Val7=/del")
except TinyHGVSError as error:
    print(error.code)
    print(error.fragment)
```
