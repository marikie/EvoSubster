# T005: Identify a behavior-preserving Stage 0 refactor

## Goal

Inspect `src/select/assembly_selection.R` and its tests, identify concrete maintainability problems, and propose the smallest high-value behavior-preserving refactor.

## Inputs and scope

- Read `.conductor/context-refactor-stage0.md` first.
- Read `src/select/assembly_selection.R` and the Stage 0 tests.
- This is read-only analysis. Do not modify production or test code.

## Constraints

- Preserve ranking order, audit/review schemas, validation errors, CLI behavior, and default audit-only QC policy.
- Prefer extraction of cohesive helpers and removal of duplicated state transitions over cosmetic renaming.
- Do not recommend a rewrite or a new dependency.
- Do not revert unrelated changes.

## Deliverable

Write `.conductor/handoffs/T005-handoff.md` using the conductor handoff contract. Include exact function/line regions, the proposed helper boundaries, regression risks, and the tests that pin behavior.

## Claim boundary

You may recommend a refactor. You must not claim it is safe or behavior-preserving until implemented and independently verified.
