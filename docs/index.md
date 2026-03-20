# tinyhgvs

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
env -u CONDA_PREFIX uv run maturin develop
```


## Quick Examples

Parse a coding DNA substitution:

```python
from tinyhgvs import parse_hgvs

variant = parse_hgvs("NM_004006.2:c.357+1G>A")
print(variant.coordinate_system.value)
print(variant.description.edit)
```

Inspect an unsupported syntax error:

```python
from tinyhgvs import TinyHGVSError, parse_hgvs

try:
    parse_hgvs("NC_000001.11:g.[123G>A;345del]")
except TinyHGVSError as error:
    print(error.code)
    print(error.fragment)
```
