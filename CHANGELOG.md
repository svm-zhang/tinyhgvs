# Changelog

All notable changes to this project will be documented in this file.

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
