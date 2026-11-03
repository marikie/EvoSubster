# Trio selection refactor implementation plan

**Goal:** Preserve selection behavior while splitting the selector by responsibility.

**Architecture:** The existing entrypoint loads seven helpers into its own environment and retains main orchestration. Existing assembly ranking and execution metadata modules remain dependencies.

**Approved specification:** User-approved plan in the 2026-10-29 conversation.

## Constraints
- Preserve CLI, function signatures, TSV types/order, selection stopping rules, caching and cleanup safety.
- Preserve existing untracked work; add no dependencies or algorithm changes.
- Only one implementation writer; parent reviews the resulting diff and evidence.

## Tasks
- [x] Extract options, IO, cleanup, Stage 0, search, cache and tree output modules.
- [x] Separate per-pair search, invalid-cache cleanup, artifact ownership, cleanup phases and annotation preparation.
- [x] Add source-environment and search-boundary regression coverage.
- [x] Run existing selector suites, deterministic before/after comparisons and hook status.

## Validation
Run test_trio_selection_helpers.R, test_outgroup_selection.R, test_tree_outputs.R,
test_selection_execution.R, test_stage0_quality_selection.R, test_stage0_gate.R,
test_trio_modules.R and test_sbst_from_dwl_cache.sh. Compare CLI help and representative
selection TSV and PDF renderings against the pre-change snapshot, without live downloads.
