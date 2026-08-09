# Stage 0 Reference-First QC Review Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make eukaryotic Stage 0 reference-first and annotation-aware, keep external QC audit-only by default, and allow only an explicit strict BUSCO override.

**Architecture:** Compute a metadata-only baseline before external QC. Attach comparable BUSCO evidence as audit and review state; only `--allow-qc-override` may replace the baseline under strict dominance rules. Preserve existing downstream whole-genome behavior when annotation is absent.

**Tech Stack:** Base R, Python 3 standard library, Bash, NCBI Datasets JSON/TSV metadata.

## Global Constraints

- Work on `feature/stage0-reference-quality-ranking`; do not touch untracked `analysis/`, `stage0-processing-flow.html`, or `stage1-matching-algorithm.html`.
- Eukaryote baseline tiers are exactly: annotated NCBI Reference; unannotated NCBI Reference; annotated RefSeq; other annotated; existing metadata fallback.
- Missing annotation is neither exclusion nor review: preserve whole-genome analysis, and use GFF-dependent non-coding analysis only when GFF exists.
- Existing prokaryote ANI/CheckM/Merqury ranking is unchanged.
- External BUSCO/Merqury never changes default selection.
- Explicit QC override may select an unannotated candidate.
- BUSCO components are percentages in `[0, 100]`; component equality tolerance is `0.2` percentage points.
- Paired GCA/GCF accessions are one assembly and cannot create a review or override.
- Use TDD: write and observe each behavior test failing before production changes.

---

### Task 1: Align the written specification and fetch paired assembly metadata

**Files:**
- Modify: `docs/superpowers/specs/2026-08-07-stage0-reference-first-qc-review-design.md`
- Modify: `src/select/fetch_assembly_metadata.py`
- Test: `test/test_fetch_metadata_taxon.py`

**Interfaces:**
- Produces metadata column `paired_accession: string`, copied from the top-level NCBI dataset report and empty when absent.
- Leaves all existing metadata columns and TSV behavior intact.

- [ ] Update the design document to the exact five-tier order in Global Constraints. Remove `no_annotated_candidate` as a review trigger and remove annotation as an override precondition.
- [ ] Add a failing Python test using a report with `paired_accession = "GCA_901000725.3"`; assert that `report_to_row()` returns that value and that `COLUMNS` still matches emitted keys.
- [ ] Run `python3 -m unittest test.test_fetch_metadata_taxon`; verify failure because `paired_accession` is absent.
- [ ] Add `paired_accession` to `COLUMNS` and `report_to_row()`, defaulting to `""`.
- [ ] Re-run the Python test and verify success.
- [ ] Run `git diff --check`, then commit only Task 1 files with `select: record paired assembly accessions`.

### Task 2: Implement the reference-first baseline and audit-only review state

**Files:**
- Modify: `src/select/assembly_selection.R`
- Test: `test/test_stage0_quality_selection.R`

**Interfaces:**
- `rank_assembly_candidates(..., allow_qc_override = FALSE)` returns `audit`, `shortlist`, `best`, and `review`.
- `select_best_assemblies(..., allow_qc_override = FALSE)` forwards the parameter.
- Audit fields: `baseline_selected`, `qc_preferred`, `review_required`, `review_reason`, `override_applied`, `assembly_equivalence_key`.

- [ ] Replace the existing eukaryote test expectation with a failing default-policy assertion: external QC does not replace the annotated Reference.
- [ ] Add failing tier tests: unannotated Reference beats annotated non-reference RefSeq; annotated RefSeq beats other annotated assembly; missing annotation alone creates no review.
- [ ] Run `Rscript test/test_stage0_quality_selection.R`; verify failures reflect the old BUSCO-first and annotation order.
- [ ] Derive `is_refseq` from `source_database == "SOURCE_DATABASE_REFSEQ"` or a `GCF_` accession, and implement the exact baseline tier order with existing metadata fields as within-tier tie-breakers.
- [ ] Compute baseline independently of external QC. Keep `selected == baseline_selected` when override is false.
- [ ] Initialize all audit fields and construct `assembly_equivalence_key` as the lexically sorted `accession|paired_accession` pair, or the accession alone when unpaired.
- [ ] Create `review` from all audit rows belonging to species marked for QC review; absence of annotation or external QC alone must not mark review.
- [ ] Re-run Stage 0 quality tests, then `Rscript test/test_stage0_gate.R`; verify success.
- [ ] Run `git diff --check`, then commit Task 2 files with `select: make Stage 0 reference-first by default`.

### Task 3: Add strict BUSCO validation, review reasons, and explicit override

**Files:**
- Modify: `src/select/assembly_selection.R`
- Test: `test/test_stage0_quality_selection.R`

**Interfaces:**
- External fields: `qc_busco_single`, `qc_busco_duplicated`, `qc_busco_internal_stop` plus existing Complete, Fragmented, Missing, Merqury fields.
- Stable review reason codes: `incomplete_busco_comparison`, `one_sided_external_qc`, `alternative_higher_single_copy`, `baseline_busco_complete_below_90`, `baseline_busco_duplicated_above_5`, `complete_gain_duplication_confounded`.
- Override selection basis: `explicit_busco_override`.

- [ ] Add failing tests for: derived Single-copy; invalid values outside 0–100; `C < D`; `abs(C - S - D) > 0.2`; `abs(C + F + M - 100) > 0.2`; legacy incomplete QC accepted but unable to override.
- [ ] Add failing Plectropomus-like test: explicit flag permits a distinct candidate with higher Single-copy and no worse D/F/M.
- [ ] Add failing Larimichthys-like test: Complete increases but Duplicated worsens, so baseline stays and `complete_gain_duplication_confounded` is recorded.
- [ ] Add failing Takifugu-like test: paired GCF/GCA accessions neither review nor override each other.
- [ ] Run the Stage 0 quality test and verify the new cases fail for missing behavior.
- [ ] Implement numeric validation and Single-copy derivation. Incomplete legacy rows remain audit-valid; validation errors identify accession and field.
- [ ] Mark `qc_preferred` by Single-copy descending, Duplicated ascending, Missing ascending, Fragmented ascending, then existing metadata order.
- [ ] Populate review reasons only for distinct assembly comparisons. Baseline warning thresholds are Complete `< 90` and Duplicated `> 5`.
- [ ] With explicit override, require same genome mode/lineage/version, complete components on both sides, distinct equivalence keys, higher Single-copy, and D/F/M each no worse. Annotation is not required.
- [ ] Re-run `Rscript test/test_stage0_quality_selection.R` and `Rscript test/test_stage0_gate.R`; verify success.
- [ ] Run `git diff --check`, then commit Task 3 files with `select: require explicit strict BUSCO overrides`.

### Task 4: Integrate CLI, review output, wrapper forwarding, and documentation

**Files:**
- Modify: `src/select/trio_selection.R`
- Modify: `src/sbst_fromDwl.sh`
- Modify: `README.md`
- Test: `test/test_trio_selection_helpers.R`
- Test: `test/test_sbst_from_dwl_cache.sh`

**Interfaces:**
- New flag `--allow-qc-override`, default false.
- Tree mode always writes `assembly_review.tsv` with headers, even when it has zero rows.

- [ ] Add failing R helper tests for help text, default false, flag parsing, and forwarding into assembly ranking.
- [ ] Add failing Bash wrapper test asserting `--allow-qc-override` is accepted and forwarded to the R selector.
- [ ] Run both tests and verify failure for missing CLI behavior.
- [ ] Implement R parsing/help, pass the flag to ranking, write `assembly_review.tsv`, and log review species count without calling unreviewed choices final.
- [ ] Implement Bash help/parsing/forwarding.
- [ ] Update README: exact baseline order; external QC audit-only by default; strict explicit override; all BUSCO components; paired accession handling; review output; BUSCO override does not prove base-level QV.
- [ ] Re-run helper and wrapper tests, then all Stage 0/Python tests from Global Constraints.
- [ ] Run `git diff --check`, then commit Task 4 files plus this plan with `docs: document audit-first Stage 0 QC workflow`.

### Task 5: Whole-change verification

**Files:** None expected.

- [ ] Run: `Rscript test/test_stage0_quality_selection.R`.
- [ ] Run: `Rscript test/test_stage0_gate.R`.
- [ ] Run: `Rscript test/test_trio_selection_helpers.R`.
- [ ] Run: `python3 -m unittest test.test_fetch_metadata_taxon`.
- [ ] Run: `python3 -m unittest test.test_trio_filter`.
- [ ] Run: `bash test/test_sbst_from_dwl_cache.sh`.
- [ ] Run: `git diff --check` and inspect `git status --short --branch`.
- [ ] Confirm unrelated untracked files are untouched and do not push.
