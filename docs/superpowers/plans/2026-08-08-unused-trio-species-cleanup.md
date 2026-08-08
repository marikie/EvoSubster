# Unused Trio Species Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove only current-run data for species excluded by tree-mode trio selection while preserving selected and pre-existing data.

**Architecture:** The downloader reports newly created paths through an optional manifest without changing stdout. The R selector owns cleanup policy because it knows both the current-run artifacts and final selected accessions; the wrapper only forwards the opt-out switch and stops when selector cleanup fails.

**Tech Stack:** Bash 3.2-compatible shell, R, TSV manifests, Python and shell integration tests.

## Global Constraints

- Work only in the isolated `feature/newick-input` worktree.
- Code and comments are English; user-facing explanations are Japanese.
- Preserve the six-field `dwl_organism.sh` stdout contract.
- Never delete a path that existed before the current run.
- Never follow a symlink or delete outside the configured root.
- Workers do not commit; the main conductor verifies and commits.

---

### Task 1: Downloader artifact ownership

**Files:**
- Modify: `src/dwl_organism.sh`
- Test: `test/test_dwl_organism.sh`

**Interfaces:**
- Produces: `--artifact-manifest PATH` with columns `accession`, `artifact_type`, `relative_path`.
- `artifact_type` is `directory_tree`, `directory`, `file`, or `symlink`.

- [x] Add failing tests for a fresh organism directory, a pre-existing directory, unchanged stdout, and empty-array retry behavior on Bash 3.2.
- [x] Run `bash test/test_dwl_organism.sh` and confirm the new assertions fail.
- [x] Implement portable before/after artifact tracking and atomic manifest output.
- [x] Run the test and `bash -n src/dwl_organism.sh` until both pass.
- [x] Write `.conductor/handoffs/T001-handoff.md` without committing.

### Task 2: Selector and wrapper cleanup policy

**Files:**
- Modify: `src/select/trio_selection.R`
- Modify: `src/sbst_fromDwl.sh`
- Modify: `README.md`
- Test: `test/test_trio_selection_helpers.R`
- Test: `test/test_sbst_from_dwl_cache.sh`
- Test: `test/test_stage0_gate.R`

**Interfaces:**
- Consumes: the Task 1 artifact manifest contract.
- Produces: `cleanup_unused_species.tsv` and `--keep-unused-species-data` in selector and wrapper tree mode.

- [x] Add failing unit tests covering selected retention, unused deletion, pre-existing retention, generated-cache cleanup, zero trios, unsafe paths, and cleanup failure.
- [x] Add failing wrapper tests for opt-out forwarding and selector-failure propagation.
- [x] Run the focused tests and confirm the new assertions fail.
- [x] Implement artifact/cache tracking, safe deletion, audit output, fatal failure propagation, and CLI help.
- [x] Fix the fetched baseline Stage 0 test so deferred `dplyr` loading is explicit.
- [x] Run focused R/shell tests until they pass.
- [x] Write `.conductor/handoffs/T002-handoff.md` without committing.

### Task 3: Adversarial verification and integration

**Files:**
- Read all Task 1 and Task 2 diffs.
- Write: `.conductor/handoffs/T003-handoff.md`

- [x] Try to construct path traversal, root deletion, symlink-following, shared-data deletion, and zero-trio counterexamples.
- [x] Check that cleanup failure prevents every downstream `sbst.sh` call.
- [x] Run focused tests and report pass/fail only.
- [x] Main conductor independently fixes verified findings and reruns all load-bearing checks.
- [x] Record the destructive-default release decision and delegation effect.
- [x] Commit documentation, implementation, tests, and the decision record locally; do not push.
