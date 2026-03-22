# Changelog

All notable changes to this project will be documented in this file.

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
