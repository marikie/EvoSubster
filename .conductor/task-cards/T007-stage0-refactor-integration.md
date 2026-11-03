# T007: Integrate and verify the Stage 0 refactor

## Goal

Apply the approved behavior-preserving extraction in `src/select/assembly_selection.R`, add characterization coverage for previously implicit override behavior, and verify the final repository state.

## Scope

- Read `.conductor/context-refactor-stage0.md` and the T005/T006 handoffs.
- Modify only `src/select/assembly_selection.R`, `test/test_stage0_quality_selection.R`, and conductor records for this phase.
- Preserve all public behavior, output columns, reason codes, validation messages, and ranking rules.
- Do not modify or stage unrelated untracked files.
- Do not push.

## Machine checks

- All Stage 0 R tests.
- Trio helper and wrapper integration tests.
- Metadata/filter Python tests.
- Shell syntax and `git diff --check`.
- Lint all new conductor handoffs.

## Claim boundary

Behavior preservation may be claimed only after the final checks pass and the independent T006 review is approved.
