# tinyhgvs

`tinyhgvs` is a lightweight HGVS parsing Python library backed by a core parser and
structured data and error model implemented in Rust.

## Installation

`tinyhgvs` can be installed directly from PyPI:

```bash
pip install tinyhgvs
```

## Quick start

- An splicing-site substitution crossing exon/intron boundary:

```python
from tinyhgvs import parse_hgvs

variant = parse_hgvs("NM_004006.2:c.357+1G>A")
print(variant.reference.primary.id)
print(variant.coordinate_system.value)
print(variant.description.location.start.coordinate)
print(variant.description.location.start.offset)
```

- Known but unsupported HGVS syntax raises `tinyhgvserror` with diagnostic code:

```python
from tinyhgvs import TinyHGVSError, parse_hgvs

try:
    parse_hgvs("NC_000001.11:g.[123G>A;345del]")
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

`tinyhgvs` is licensed under MIT license.
