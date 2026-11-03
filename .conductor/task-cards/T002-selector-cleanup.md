# Task card: T002 - Selector cleanup and wrapper propagation

## Goal

Implement current-run artifact/cache cleanup after trio selection, plus the
opt-out flag and failure propagation, using TDD.

## Authority / inputs

- Context: `.conductor/context-packet.md`
- Design: `docs/superpowers/specs/2026-08-08-unused-trio-species-cleanup-design.md`
- Manifest contract from T001 is fixed by the context packet.
- You own `src/select/trio_selection.R`, `src/sbst_fromDwl.sh`, `README.md`,
  `test/test_trio_selection_helpers.R`, `test/test_sbst_from_dwl_cache.sh`, and
  `test/test_stage0_gate.R`.

## CLAIM BOUNDARY (hard)

- Allowed changes: selector-owned cleanup, wrapper opt-out forwarding, focused
  documentation/tests, and the existing deferred-dplyr baseline test regression.
- MUST NOT edit downloader files or unrelated pipeline code.
- No git commit. Do not revert other workers' changes.

## Steps

1. Add failing cleanup helper and wrapper integration tests.
2. Verify RED, implement minimal behavior, verify GREEN.
3. Ensure selected and pre-existing artifacts survive, zero-trio cleanup works,
   unsafe paths fail, and selector failure blocks downstream execution.
4. Update CLI help/README and write the compact handoff.

## Machine checks

- `Rscript test/test_trio_selection_helpers.R`
- `Rscript test/test_stage0_gate.R`
- `Rscript test/test_outgroup_selection.R`
- `bash test/test_sbst_from_dwl_cache.sh`
- `bash -n src/sbst_fromDwl.sh`
- `git diff --check -- src/select/trio_selection.R src/sbst_fromDwl.sh README.md test/test_trio_selection_helpers.R test/test_sbst_from_dwl_cache.sh test/test_stage0_gate.R`

## Stop conditions

- Stop on unsafe deletion ambiguity or any failing focused check.

## Handoff

Write `.conductor/handoffs/T002-handoff.md` with no raw logs.
