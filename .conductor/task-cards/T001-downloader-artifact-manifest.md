# Task card: T001 - Downloader artifact manifest

## Goal

Implement the optional artifact manifest in the downloader using TDD, including
the fetched Bash 3.2 empty-array regression.

## Authority / inputs

- Context: `.conductor/context-packet.md`
- Design: `docs/superpowers/specs/2026-08-08-unused-trio-species-cleanup-design.md`
- You own only `src/dwl_organism.sh` and `test/test_dwl_organism.sh`.

## CLAIM BOUNDARY (hard)

- Allowed changes: downloader manifest behavior and its tests.
- MUST NOT change selector, wrapper, README, or other tests.
- Preserve the downloader's six-field stdout contract.
- No git commit. Do not revert other workers' changes.

## Steps

1. Add failing assertions for fresh and pre-existing output directories.
2. Verify RED, implement minimal Bash 3.2-compatible behavior, verify GREEN.
3. Write an atomic TSV manifest with only newly created relative paths.
4. Run all machine checks and write the compact handoff.

## Machine checks

- `bash test/test_dwl_organism.sh`
- `bash -n src/dwl_organism.sh`
- `git diff --check -- src/dwl_organism.sh test/test_dwl_organism.sh`

## Stop conditions

- Stop on any failing check or contract ambiguity.

## Handoff

Write `.conductor/handoffs/T001-handoff.md` with no raw logs.
