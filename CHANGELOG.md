# Changelog

All notable changes to this project will be documented in this file.

## [0.6.2]

### Uncertain location support
- Added support for uncertain location syntax across DNA, RNA, and protein
  descriptions, including:
  - DNA uncertain single-region locations such as
    `NC_000023.10:g.(33038277_33038278)C>T`
  - DNA variant with uncertain left and right breakpoints such as
    `NC_000023.10:g.(?_234567)_(345678_?)del`
  - RNA uncertain locations such as `NM_004006.2:r.(71_72)_(90_91)del`
  - Protein uncertain locations such as `p.(Ala123_Pro131)Ter`

### Rust and Python models
- Added a shared Rust `Location<T>` model for known and uncertain locations.
- Updated nucleotide and protein descriptions to use `Location<T>` for edited
  locations.
- Added Rust helpers for known starts/ends and uncertain left/right regions.
- Mirrored the location model through the PyO3 bridge and Python surface.

### Diagnostics and tests
- Moved uncertain location syntax out of the unsupported boundary.
- Kept uncertain size syntax under `unsupported.uncertain_size`.
- Kept uncertain protein consequence forms under
  `unsupported.protein_uncertain_consequence`.
- Added Rust and Python coverage for valid DNA, RNA, and protein uncertain
  locations, malformed uncertain locations, and remaining unsupported uncertain
  syntax.

### Documentation and support inventory
- Updated the unsupported syntax inventory to mark uncertain location support
  as available since `0.6.2`.
- Updated Python docstrings for the location model and coordinate state helpers.

## [0.6.1]

### Protein allele support
- Added support for protein allele syntax on the current HGVS path, including:
  - one written allele such as `p.[Ser68Arg;Asn594del]`
  - two written alleles in trans such as `p.[Ser68Arg];[Ser68=]`
  - uncertain phase such as `p.(Ser73Arg)(;)(Asn103del)`
  - Predicted protein allele form such as `p.[(Ser73Arg;Asn103del)]`
  - Allele with mixed known and predicted variants such as `p.[Phe233Leu;(Cys690Trp)]`
  - `p.[Ser86Arg];[0]`
- Retired the `unsupported.protein_allele` family and replaced with more grained
   code for protein allele syntax that are yet to be supported (see below).

### Rust and Python models
- Reused the shared `AlleleVariant`, `Allele`, and `AllelePhase` model shape
  for protein allele parsing on both the Rust and Python sides.
- Added `VariantDescription::ProteinAllele(...)` on the Rust side and mirrored
  the same shape through the PyO3 bridge and Python surface.
- Moved the shared Python allele model into the shared module.

### Diagnostics and tests
- Kept the remaining out-of-scope protein allele forms under structured
  diagnostics, including:
  - `unsupported.allele_unknown_variant`
  - `unsupported.alternate_allele_state`
  - `unsupported.one_allele_multi_protein`
- Added Rust and Python coverage for both valid and malformed protein allele
  examples.

### Documentation and support inventory
- Updated the unsupported syntax inventory to mark protein allele support as
  available since `0.6.1`.
- Refreshed shared Python allele docstrings and public error examples so they
  reflect the current protein allele boundary.

## [0.6.0]

### Nucleotide allele support
- Added support for supported DNA and RNA allele syntax, including:
  - single-allele cis forms such as `g.[123G>A;345del]`
  - variants *in trans* such as `r.[123c>a];[345del]`
  - uncertain-phase forms such as `g.123G>A(;)345del`
  - mixed-phase chains such as `c.[296T>G;476T>C];[476T>C](;)1083A>C`
- Removed supported DNA/RNA allele syntax from the unsupported parser boundary.

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
- Added support for supported HGVS protein frameshift syntax, including short form
  such as `p.Arg97fs` and long form such as `p.Arg97ProfsTer23` and `p.Arg97ProfsTer?`.
- Removed supported protein frameshift syntax from the unsupported parser boundary
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
- Removed supported repeat families from the unsupported parser boundary while
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
