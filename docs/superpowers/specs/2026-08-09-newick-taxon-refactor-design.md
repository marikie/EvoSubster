# Newick Taxon Processing Refactor Design

## Goal

Refactor the accession-free Newick conversion and taxon-name normalization code
without changing any observable behavior, command-line interface, output schema,
message text, or conversion result.

## Scope

The refactor is limited to:

- `src/select/strip_newick_accessions.py`
- the taxon-label normalization helper in `src/select/trio_selection.R`
- the organism-name normalization helper in
  `src/select/fetch_assembly_metadata.py`
- directly related tests

Stage 0 ranking, QC policy, tree pruning, Stage 1 matching, wrapper options,
README content, and output TSV schemas are out of scope.

## Design

### Newick converter

Keep the converter in one module, but separate its responsibilities clearly:

1. `scan_leaf_segment(text, start)` consumes comments, quoted labels, and one
   complete leaf segment.
2. An immutable `LeafScan` result carries the segment end, logical label,
   source-character positions, and whether a label was present.
3. A rewrite function removes only the source positions corresponding to a
   terminal versioned GCA/GCF accession.
4. An immutable `ConversionResult` carries the rewritten text, converted count,
   and total leaf count. `rewrite_newick(text)` returns this result.
5. `convert_file(input_path, output_path, overwrite)` owns input validation,
   reading, atomic writing, and returns the conversion result to the CLI adapter.
6. `main()` remains responsible only for argument parsing, stable diagnostics,
   and exit codes.

The current two-pass `consume_leaf_segment()` plus `rewrite_leaf_segment()`
logic will become one leaf scan. Comments and Newick quoting remain byte-for-byte
preserved, while only accession characters are removed.

Atomic no-clobber behavior remains unchanged: without `--force`, publishing the
completed temporary file uses an exclusive hard link; with `--force`, the
temporary file atomically replaces the destination.

### Taxon normalization

In R, extract the scalar normalization rules into `taxon_from_label()` and keep
`species_from_label()` as the vectorized adapter. The rules remain unchanged:
optional Newick outer quotes are decoded, actual whitespace is rejected,
versioned or versionless terminal GCA/GCF accessions are rejected, all
underscore-delimited components are required and retained, and at least two
components are required.

In Python, extract whitespace normalization into `normalize_organism_name()` and
keep `split_organism_name()` as the dictionary-producing adapter. `species`
remains the complete normalized NCBI organism name and `genus` remains its first
token.

## Compatibility Invariants

- The converter CLI remains `--input`, `--output`, and optional `--force`.
- Exit codes and all current user-facing messages remain identical.
- Output Newick text remains byte-identical to the current implementation.
- Comments, quoting, internal-node labels, whitespace, topology, and branch
  lengths remain unchanged.
- The `species` and `genus` values produced by R and Python remain identical to
  the current implementation for every accepted and rejected test input.
- No TSV columns or wrapper options change.

## Verification

Existing characterization tests remain the primary contract. Add focused unit
tests only where the new named result types or extracted scalar helpers need a
direct contract.

Behavior preservation will also be checked by comparing the refactored code
against commit `bbcaad8` over deterministic randomized Newick strings covering
plain and quoted leaves, escaped quotes, comments before/inside/after labels,
internal-node labels, mixed converted/unconverted leaves, whitespace, branch
lengths, and both GCA/GCF suffixes.

The final verification includes the full Python suite, all Stage 0 and outgroup
R suites, both shell suites, Python and shell syntax checks, `git diff --check`,
and the Eupercaria conversion/Stage 0 dry-run used for the original feature.

## Git History

The design is committed separately. The behavior-preserving implementation and
its test adjustments will be a second focused commit. Existing untracked HTML,
analysis output, Conductor files, and the preceding `.gitignore` commit remain
untouched.
