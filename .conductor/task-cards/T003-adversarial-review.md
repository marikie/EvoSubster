# Task card: T003 - Adversarial cleanup review

## Goal

Read the integrated diff and try to refute its safety and completeness claims.

## Authority / inputs

- Context: `.conductor/context-packet.md`
- Design and implementation diff in the isolated worktree.
- Read-only review; do not edit source or tests.

## CLAIM BOUNDARY (hard)

- Allowed: report concrete findings with paths/lines and machine-check results.
- MUST NOT modify files, commit, or promote the change as release-ready.

## Steps

1. Test path traversal, root deletion, symlink handling, and shared-data cases.
2. Verify selected, unused, zero-trio, cache, opt-out, and failure propagation.
3. Run focused tests and write a compact handoff.

## Machine checks

- `bash test/test_dwl_organism.sh`
- `Rscript test/test_trio_selection_helpers.R`
- `bash test/test_sbst_from_dwl_cache.sh`
- `git diff --check`

## Stop conditions

- Escalate any ambiguous destructive behavior; do not assume safety.

## Handoff

Write `.conductor/handoffs/T003-handoff.md` with no raw logs.
