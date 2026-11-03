# Newick Taxon Processing Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Simplify the Newick converter and R/Python taxon normalization helpers without changing observable behavior.

**Architecture:** Replace the converter's two leaf-segment passes and boolean tuples with one scanner and immutable named results. Extract scalar normalization helpers in R and Python while retaining the current public adapters, CLI, diagnostics, Newick output, TSV schema, and Stage 0 join behavior.

**Tech Stack:** Python 3 standard library, base R, Bash, `unittest`, `ape` for end-to-end verification.

## Global Constraints

- Work on `feature/stage0-reference-quality-ranking` and preserve the existing `.gitignore` commit.
- Do not touch untracked HTML, `analysis/`, or Conductor files.
- Do not change converter arguments, exit codes, diagnostic text, output bytes, or atomic no-clobber behavior.
- Do not change accepted/rejected taxon labels, NCBI exact-match behavior, TSV columns, Stage 0 ranking, tree pruning, or Stage 1 matching.
- Keep the converter in `src/select/strip_newick_accessions.py`; do not create a parser package.
- Use characterization tests plus TDD for every new interface.
- Produce one focused implementation commit after all refactor tasks and verification pass, as required by the approved design.

---

### Task 1: Refactor the Newick scanner and file-conversion boundary

**Files:**
- Modify: `src/select/strip_newick_accessions.py`
- Test: `test/test_strip_newick_accessions.py`

**Interfaces:**
- `LeafScan(end: int, logical_label: str, source_positions: Tuple[int, ...], has_label: bool)` is a frozen dataclass.
- `ConversionResult(text: str, converted: int, leaf_count: int)` is a frozen dataclass.
- `CliUsageError(ValueError)` identifies current exit-code-2 path errors.
- `scan_leaf_segment(text: str, start: int) -> LeafScan` replaces `consume_leaf_segment()` and the scanning half of `rewrite_leaf_segment()`.
- `rewrite_leaf_segment(text: str, start: int, scan: LeafScan) -> Tuple[str, bool]` removes only matched source positions.
- `rewrite_newick(text: str) -> ConversionResult` replaces the three-item tuple.
- `convert_file(input_path: Path, output_path: Path, overwrite: bool) -> ConversionResult` validates paths, reads, rewrites, and publishes atomically.

- [ ] **Step 1: Add failing tests for the named scanner and conversion interfaces**

Add direct tests equivalent to:

```python
def test_scan_leaf_segment_exposes_logical_label_and_source_positions(self):
    text = "A_a[inside]_GCA_000000001.1:0.2"
    scan = CONVERTER.scan_leaf_segment(text, 0)
    self.assertEqual(scan.end, text.index(":"))
    self.assertEqual(scan.logical_label, "A_a_GCA_000000001.1")
    self.assertTrue(scan.has_label)
    self.assertEqual(
        "".join(text[position] for position in scan.source_positions),
        scan.logical_label,
    )

def test_rewrite_newick_returns_named_counts(self):
    result = CONVERTER.rewrite_newick("(A_a_GCA_000000001.1:1,B_b:2);\n")
    self.assertIsInstance(result, CONVERTER.ConversionResult)
    self.assertEqual(result.text, "(A_a:1,B_b:2);\n")
    self.assertEqual((result.converted, result.leaf_count), (1, 2))

def test_convert_file_owns_path_validation(self):
    with self.assertRaisesRegex(
        CONVERTER.CliUsageError,
        "--input and --output must be different paths",
    ):
        CONVERTER.convert_file(tree_path, tree_path, overwrite=True)
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
python3 -m unittest test.test_strip_newick_accessions
```

Expected: failures because `LeafScan`, `ConversionResult`, `CliUsageError`,
`scan_leaf_segment()`, and `convert_file()` do not exist and
`rewrite_newick()` still returns a tuple.

- [ ] **Step 3: Implement the immutable result types and one-pass leaf scanner**

Add:

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class LeafScan:
    end: int
    logical_label: str
    source_positions: Tuple[int, ...]
    has_label: bool

@dataclass(frozen=True)
class ConversionResult:
    text: str
    converted: int
    leaf_count: int

class CliUsageError(ValueError):
    """Invalid CLI path configuration that retains exit code 2."""
```

Implement `scan_leaf_segment()` as one cursor pass. It must skip comments from
the logical label, decode doubled quotes as one logical apostrophe, retain
absolute source positions for every logical character, stop only at an
unquoted/uncommented `(`, `)`, `:`, `,`, or `;`, and preserve unterminated
comment/quote errors from the existing consumers.

Implement `rewrite_leaf_segment()` by matching `ACCESSION_SUFFIX_RE` against
`scan.logical_label` and removing the corresponding absolute positions from
`text[start:scan.end]`. Keep comments, whitespace, and quote characters.

- [ ] **Step 4: Make `rewrite_newick()` and the CLI use named results**

Return:

```python
return ConversionResult("".join(output), converted, leaf_count)
```

Implement `convert_file()` so same-path, missing-input, missing-output-directory,
and destination-exists cases raise `CliUsageError` with the exact existing text.
Keep `write_text_atomic()` unchanged. Map `CliUsageError` to exit code 2 in
`main()` and all other current conversion errors to exit code 1. Keep the two
success diagnostics byte-identical.

- [ ] **Step 5: Run the converter tests and verify GREEN**

Run:

```bash
python3 -m unittest test.test_strip_newick_accessions
python3 -m py_compile src/select/strip_newick_accessions.py
```

Expected: all converter tests pass and compilation exits 0.

---

### Task 2: Extract scalar taxon normalization helpers

**Files:**
- Modify: `src/select/trio_selection.R`
- Modify: `src/select/fetch_assembly_metadata.py`
- Test: `test/test_trio_selection_helpers.R`
- Test: `test/test_fetch_metadata_taxon.py`

**Interfaces:**
- `taxon_from_label(name) -> character(1)` owns the current scalar R normalization rules.
- `species_from_label(labels) -> character(length(labels))` remains the vector adapter with `USE.NAMES = FALSE`.
- `normalize_organism_name(organism_name: str) -> str` collapses whitespace only.
- `split_organism_name(organism_name: str) -> Dict[str, str]` remains the dictionary adapter.

- [ ] **Step 1: Add failing direct helper tests**

Add an R assertion equivalent to:

```r
check(
  identical(taxon_from_label("'Chaunax_sp._Z400'"), "Chaunax sp. Z400") &&
    is.na(taxon_from_label("Chaunax_sp._Z400_GCA_037577475")),
  "scalar taxon parser preserves accepted and rejected label behavior"
)
```

Add a Python assertion equivalent to:

```python
def test_normalizes_complete_organism_name_whitespace(self):
    self.assertEqual(
        f.normalize_organism_name("  Genus\talpha  strain-X "),
        "Genus alpha strain-X",
    )
```

- [ ] **Step 2: Run helper tests and verify RED**

Run:

```bash
Rscript test/test_trio_selection_helpers.R
python3 -m unittest test.test_fetch_metadata_taxon
```

Expected: failures because the two scalar helper names do not exist.

- [ ] **Step 3: Extract helpers without changing adapters**

Move the current scalar R body into `taxon_from_label()` and implement:

```r
species_from_label <- function(labels) {
  vapply(labels, taxon_from_label, character(1), USE.NAMES = FALSE)
}
```

Implement Python normalization as:

```python
def normalize_organism_name(organism_name: str) -> str:
    return " ".join(organism_name.split())

def split_organism_name(organism_name: str) -> Dict[str, str]:
    species = normalize_organism_name(organism_name)
    tokens = species.split()
    genus = tokens[0] if tokens else ""
    return {"genus": genus, "species": species}
```

- [ ] **Step 4: Run helper and Stage 0 tests and verify GREEN**

Run:

```bash
Rscript test/test_trio_selection_helpers.R
Rscript test/test_stage0_gate.R
Rscript test/test_stage0_quality_selection.R
python3 -m unittest test.test_fetch_metadata_taxon
```

Expected: all tests pass with unchanged Stage 0 results.

---

### Task 3: Prove behavior preservation, review, and commit

**Files:**
- Modify only if review finds a defect: files from Tasks 1-2 and their tests

**Interfaces:** No additional interfaces.

- [ ] **Step 1: Run deterministic baseline/current differential comparison**

Load `src/select/strip_newick_accessions.py` from commit `bbcaad8` into a
temporary module and compare it with the worktree module for at least 2,000
seeded cases. Normalize the current result to
`(result.text, result.converted, result.leaf_count)` before comparison. Generate
trees from combinations of:

```text
plain leaves
single-quoted leaves
doubled apostrophes inside quoted leaves
comments before, inside, and after leaf-label components
quoted and unquoted internal-node labels
GCA and GCF versioned suffixes
versionless and accession-free leaves
spaces and newlines around structural delimiters
integer, decimal, and scientific-notation branch lengths
```

Use random seed `20260809`. Expected: 2,000/2,000 rewritten text and count
triples are identical.

- [ ] **Step 2: Run the Eupercaria structural and Stage 0 dry-run check**

Convert
`/Users/marikonakagawa/biohazard/TheTreeOfLife/Eupercaria_tree.nwk` in a
temporary directory. Verify:

```text
Converted 193 of 193 leaf labels
0 terminal versioned GCA/GCF suffixes remain
Ntip = 193
edge, edge.length, Nnode, and node.label are unchanged
tree_taxa.txt contains Chaunax sp. Z400
assembly_metadata.tsv contains Chaunax sp. Z400
selected_assemblies.tsv selects GCA_037577475.1 for Chaunax sp. Z400
```

- [ ] **Step 3: Run the complete regression suite**

Run:

```bash
git diff --check
python3 -m py_compile \
  src/select/strip_newick_accessions.py \
  src/select/fetch_assembly_metadata.py
bash -n src/sbst_fromDwl.sh
(cd test && python3 -m unittest discover -p 'test_*.py')
Rscript test/test_trio_selection_helpers.R
Rscript test/test_stage0_gate.R
Rscript test/test_stage0_quality_selection.R
Rscript test/test_outgroup_selection.R
bash test/test_sbst_from_dwl_cache.sh
bash test/test_train_cache.sh
```

Expected: every command exits 0. The pre-existing `test_subRatio.py`
unclosed-file `ResourceWarning` may still appear but must not create a test
failure.

- [ ] **Step 4: Request an independent behavior-preservation review**

Provide the reviewer with commit `bbcaad8`, the approved design, the complete
working-tree diff, differential-test evidence, and Eupercaria evidence. Fix all
Critical, Important, and Minor findings, then rerun the affected tests.

- [ ] **Step 5: Commit only the refactor implementation**

Run:

```bash
git add \
  src/select/strip_newick_accessions.py \
  src/select/trio_selection.R \
  src/select/fetch_assembly_metadata.py \
  test/test_strip_newick_accessions.py \
  test/test_trio_selection_helpers.R \
  test/test_fetch_metadata_taxon.py
git diff --cached --check
git commit -m "refactor: simplify Newick taxon processing"
```

Confirm the implementation commit does not include `.gitignore`, untracked
HTML, `analysis/`, or Conductor files. Do not push unless explicitly requested.
