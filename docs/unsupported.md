# Unsupported Syntax

`tinyhgvs` recognizes a number of HGVS syntax families that are not yet modeled
for successful parsing in the current release. When one of these families is
detected, the parser raises a structured error:

- Python: `TinyHGVSError.code`
- Rust: `ParseHgvsError.code`

All unsupported diagnostic codes exposed by tinyhgvs begin with the prefix `unsupported.`.
For readability, the tables below show only the category suffix.

The tables below track syntax families by molecule type. Rows with `-` in the
`Supported since` column are still unsupported in the current release. Rows
with a version number are retained to show when a previously cataloged syntax
family became supported.

Each row shows one representative example. Additional examples are grouped in
the collapsible sections below each table.

## DNA

| Diagnostic code category | Biology category | Representative example | Supported since |
| --- | --- | --- | --- |
| `allele` | DNA allele | `NC_000001.11:g.[123G>A;345del]` | `0.6.0` |
| `allele_unknown_variant` | DNA allele with variant unknown | `NC_000001.11:g.[123G>A];[?]` | `-` |
| `allele_uncertain_variant_state` | DNA allele with uncertain variant state | `NC_000001.11:g.[123G>A](;)(345del)` | `-` |
| `uncertain_range` | DNA range-uncertain edit | `NC_000023.10:g.(33038277_33038278)C>T` | `-` |
| `dna_repeat` | DNA repeated sequence | `NC_000014.8:g.123CAG[23]` | `0.2.0` |
| `telomeric_position` | DNA telomeric coordinate | `NC_000023.11:g.pter_qtersup` | `-` |
| `epigenetic_edit` | DNA epigenetic edit | `NC_000011.10:g.1999904_1999946|gom` | `-` |
| `cdna_offset_anchor` | UTR-anchored intron edit | `NM_001385026.1:c.-666+629C>T` | `0.5.0` |


??? Example "More examples"

    === "allele"

        - `NC_000001.11:g.[123G>A];[345del]`
        - `NC_000001.11:g.123G>A(;)345del`

    === "allele_unknown_variant"

        - `NC_000001.11:g.[123G>A;345del];[?]`

    === "allele_uncertain_variant_state"

        - `NC_000001.11:g.[123G>A](;)(345del)`

    === "uncertain_range"

        - `NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup`
        - `NC_000023.11:g.(?_31120496)_(33339477_?)del`

    === "telomeric_position"

        - `NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]`

## RNA

| Diagnostic code category | Biology category | Representative example | Supported since |
| --- | --- | --- | --- |
| `allele` | RNA allele | `LRG_199t1:r.[76a>u;103del]` | `0.6.0` |
| `allele_unknown_variant` | RNA allele with variant unknown | `NM_004006.3:r.[123c>a];[?]` | `-` |
| `allele_uncertain_variant_state` | RNA allele with uncertain variant state | `NM_004006.3:r.[123c>a](;)(345del)` | `-` |
| `rna_special_state` | RNA special-state outcome | `NM_004006.3:r.?` | `-` |
| `rna_uncertain_position` | RNA uncertain insertion | `NM_004006.2:r.(222_226)insg` | `-` |
| `rna_repeat` | RNA repeated sequence | `NM_004006.3:r.-124_-123[14]` | `0.2.0` |
| `rna_splicing_outcome` | RNA uncertain splicing | `NC_000023.11(NM_004006.2):r.spl` | `-` |
| `rna_adjoined_transcript` | RNA adjoined transcript | `NM_002354.2:r.-358_555::NM_000251.2:r.212_*279` | `-` |


??? Example "More examples"

    === "allele"

        - `NM_004006.3:r.[123c>a;345del]`
        - `NM_004006.3:r.123c>a(;)345del`

    === "allele_unknown_variant"

        - `NM_004006.3:r.[123c>a;345del];[?]`

    === "allele_uncertain_variant_state"

        - `NM_004006.3:r.[123c>a](;)(345del)`

    === "rna_special_state"

        - `NM_004006.3:r.(1388g>a)`
        - `NM_004006.3:r.0`
        - `NM_004006.3:r.spl`

    === "rna_splicing_outcome"

        - `NC_000023.11(NM_004006.2):r.[897u>g,832_960del]`
        - `NC_000023.11(NM_004006.2):r.?`

    === "rna_adjoined_transcript"

        - `NM_152263.2:r.-115_775::aggcucccuugg::NM_002609.3:r.1580_*1924`

## Protein


| Diagnostic code category | Biology category | Representative example | Supported since |
| --- | --- | --- | --- |
| `allele` | Protein allele | `NP_003997.1:p.Val7=/del` | `-` |
| `protein_frameshift` | Protein frameshift | `NP_0123456.1:p.Arg97fs` | `0.3.0` |
| `protein_extension` | Protein extension | `NP_003997.2:p.Met1ext-5` | `0.4.0` |
| `protein_repeat` | Protein repeated sequence | `NP_0123456.1:p.Ala2[10]` | `0.2.0` |
| `protein_insertion_payload` | Unknown/truncating protein insertion | `p.Arg78_Gly79insXaa[23]` | `-` |
| `protein_uncertain_consequence` | Protein uncertain consequence | `p.(Gly719Ala^Ser)` | `-` |


??? Example "More examples"

    === "allele"

        - `LRG_199p1:p.Trp24=/Cys`
        - `NP_003997.1:p.[(Ser73Arg;Asn103del)]`

    === "protein_insertion_payload"

        - `NP_060250.2:p.Gln746_Lys747ins*63`
        - `NP_003997.1:p.(Val582_Asn583insXaa[5])`

    === "protein_uncertain_consequence"

        - `NP_003997.1:p.(Gly56Ala^Ser^Cys)`
