# Changelog

All notable changes to this project will be documented in this file.

## [0.3.3]

This patch addresses the security alerts from dependabot.

### Dependency maintenance
- Updated the PyO3 dependency used by the Python extension to a patched
  release line.
- Refreshed the Python docs/dev lockfile so the flagged `Pygments` dependency
  resolves to a patched version.
- Updated the docs tooling dependency resolution so `zensical build` works
  cleanly with the patched `Pygments` release line.

## [0.3.2]

### Python bridge cleanup
- Removed the unused Python-to-Rust reverse conversion path from the PyO3
  extension layer.
- Removed the internal `_roundtrip_variant` helper and the corresponding stub
  entry from the Python package.

## [0.3.1]

### Test suite cleanup
- Reworked Rust and Python parser tests so they start from raw HGVS strings and
  assert on parser-produced objects instead of comparing against hand-built
  expected model objects.
- Removed low-value tests that focused on helper or bridge behavior.

### Python test tooling
- Added persistent `pytest` coverage settings in `py-tinyhgvs/pyproject.toml`.
- Enforced a minimum Python test coverage threshold of 95%.


## [0.3.0]

### Protein frameshift support
- Added support for exact HGVS protein frameshift syntax, including short form
  such as `p.Arg97fs` and long form such as `p.Arg97ProfsTer23` and `p.Arg97ProfsTer?`.
- Removed exact protein frameshift syntax from the unsupported parser boundary
  while keeping adjacent unsupported protein families unchanged.

### Rust and Python models
- Added first-class protein frameshift models on both the Rust and Python sides.
- Extended the PyO3 bridge so frameshift variants are mapped cleanly between
  the Rust parser result and the public Python API.
- Added parser-facing valid and malformed frameshift tests.

### Documentation and support inventory
- Updated the unsupported syntax inventory to mark protein frameshift as
  supported since `0.3.0`.
- Refreshed Rust docs, Python docstrings, README examples, and rendered docs
  to present protein frameshift as a supported syntax family.

## [0.2.0]

### Repeat support
- Added support for primary HGVS repeat syntax across DNA, RNA, and protein
  variants.
- Added top-level repeat edit models on the Rust and Python sides, including
  nucleotide repeat blocks and protein repeat counts.
- Removed exact repeat families from the unsupported parser boundary while
  keeping repeat-adjacent allele and uncertain-count diagnostics explicit.

### Python and typing
- Added bidirectional PyO3 conversion support for repeat variants.
- Added public Python roundtrip coverage for repeat-bearing variants.
- Added `mypy` for the public Python package and shipped `py.typed`.

### Documentation and support inventory
- Updated the unsupported syntax inventory to mark repeat families as supported
  since `0.2.0`.
- Refreshed public examples and release-facing documentation for repeat
  support.

## [0.1.0]

### Rust
- Initial release of the `tinyhgvs` crate for parsing HGVS expressions into
    structured data.
- Published the crate on crates.io and Rust API documentation on docs.rs.

### Python
- Initial release of the `tinyhgvs` Python package backed by the Rust parser
    via PyO3.
- Published Python API documentation on the project documentation site.

### Core functionality
- Added typed HGVS data models for parsed variants.
- Added structured diagnostics that distinguish invalid syntax from
    recognized but unsupported HGVS families.
