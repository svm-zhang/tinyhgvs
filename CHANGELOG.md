# Changelog

All notable changes to this project will be documented in this file.

## [0.6.0]

### Nucleotide allele support
- Added support for exact DNA and RNA allele syntax, including:
  - single-allele cis forms such as `g.[123G>A;345del]`
  - variants *in trans* such as `r.[123c>a];[345del]`
  - uncertain-phase forms such as `g.123G>A(;)345del`
  - mixed-phase chains such as `c.[296T>G;476T>C];[476T>C](;)1083A>C`
- Removed exact DNA/RNA allele syntax from the unsupported parser boundary.

### Rust and Python models
- Added allele container model on Rust backend and mirrored in Python surface:
  - `AlleleVariant`
  - `Allele`
  - `AllelePhase`
- Exposed helper views of alleles, phase, variants written in the description.

### Diagnostics and tests
- Added specific allele diagnostics for unsupported cases that still remain out
  of scope, including `[?]` and `(;)(...)` forms.
- Added test coverage for parsing both valid and malformed allele on both Rust
  Python sides.

### Documentation and support inventory
- Updated the unsupported syntax inventory to mark DNA and RNA allele support
  as available since `0.6.0`.
- Refreshed public docs, Python docstrings, and error examples so they reflect
  the supported nucleotide allele surface.

## [0.5.0]

### Coding-DNA coordinate support
- Added support for CDS-anchored intronic coding-DNA coordinates such as
  `c.-106+2`, `c.-666+629`, `c.*639-1`, and `c.*24-12888`.
- Added support for intervals built from CDS-anchored intronic positions, such
  as `c.-490-342_-490-341del`.
- Removed `unsupported.cdna_offset_anchor` from the unsupported parser
  boundary.

### Coordinate helper semantics
- Refined nucleotide coordinate helper semantics so `is_intronic` follows any
  nonzero offset, while `is_five_prime_utr` and `is_three_prime_utr` remain
  limited to exonic UTR positions.
- Added explicit CDS-anchor helpers for nucleotide coordinates on the Rust and
  Python public models.

### Tests and support inventory
- Added test coverage on representative valid and malformed CDS-anchored
  coordinate forms on both Rust and Python side.
- Updated the unsupported syntax inventory to mark `cdna_offset_anchor` as
  supported since `0.5.0`.

## [0.4.0]

### Protein extension support
- Added support for HGVS protein extension syntax, including N-terminal
  extension such as `p.Met1ext-5` and C-terminal extension such as
  `p.Ter110GlnextTer17` and `p.Ter327ArgextTer?`.
- Accepted both `Ter` and `*` stop-codon spellings in protein extension input
  while normalizing the parsed public model to `Ter`.
- Removed `unsupported.protein_extension` code from the unsupported inventory.

### Rust and Python models
- Added protein extension models on both the Rust and Python sides.
- Extended the PyO3 bridge so protein extension variants map cleanly from the
  Rust parser result into the public Python API.
- Added test coverage for both valid and malformed extension variants.

### Documentation and support inventory
- Updated the unsupported syntax inventory to mark protein extension as
  supported since `0.4.0`.
- Refreshed Rust docs, Python docstrings, and rendered docs to present
  protein extension as a supported syntax family.

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
