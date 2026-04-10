# tinyhgvs

[![PyPI](https://img.shields.io/pypi/v/tinyhgvs.svg)](https://pypi.org/project/tinyhgvs/)
[![Python Versions](https://img.shields.io/pypi/pyversions/tinyhgvs.svg)](https://pypi.org/project/tinyhgvs/)
[![CI](https://github.com/svm-zhang/tinyhgvs/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/svm-zhang/tinyhgvs/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/svm-zhang/tinyhgvs.svg)](https://github.com/svm-zhang/tinyhgvs/blob/main/LICENSE)
[![crates.io](https://img.shields.io/crates/v/tinyhgvs.svg)](https://crates.io/crates/tinyhgvs)
[![docs.rs](https://img.shields.io/docsrs/tinyhgvs)](https://docs.rs/tinyhgvs)

`tinyhgvs` is a lightweight HGVS parsing Python library backed by a core parser and
structured data and error model implemented in Rust.

## Installation

`tinyhgvs` can be installed directly from PyPI:

```bash
pip install tinyhgvs
```

## Quick start

- A splicing-site substitution crossing an exon/intron boundary:

```python
from tinyhgvs import parse_hgvs

variant = parse_hgvs("NM_004006.2:c.357+1G>A")
print(variant.reference.primary.id)
print(variant.coordinate_system.value)
print(variant.description.location.start.coordinate)
print(variant.description.location.start.offset)
```

- An exact repeat is parsed into a structured repeat edit:

```python
from tinyhgvs import parse_hgvs

variant = parse_hgvs("NC_000014.8:g.123CAG[23]")
print(variant.description.location.start.coordinate)
print(variant.description.edit.blocks[0].unit)
print(variant.description.edit.blocks[0].count)
```

- A protein frameshift variant:

```python
from tinyhgvs import parse_hgvs

variant = parse_hgvs("NP_0123456.1:p.Arg97ProfsTer23")
print(variant.description.effect.edit.to_residue)
print(variant.description.effect.edit.stop.ordinal)
```

- Known but unsupported HGVS syntax raises `TinyHGVSError` with a diagnostic code:

```python
from tinyhgvs import TinyHGVSError, parse_hgvs

try:
    parse_hgvs("NP_003997.1:p.[Lys31Asn,Val25_Lys31del]")
except TinyHGVSError as error:
    print(error.code)
```

## Documentation

Please refer to official documents for details about both Python and Rust APIs:

- [Python API and project docs](https://svm-zhang.github.io/tinyhgvs)
- [Rust API docs](https://docs.rs/tinyhgvs)

## Citation

If you use `tinyhgvs` in your project, please kindly cite the github repository.

## License

`tinyhgvs` is licensed under the MIT license.
