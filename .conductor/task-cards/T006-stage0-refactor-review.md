# T006: Adversarial review of the Stage 0 refactor

## Goal

After implementation is available, independently inspect the final diff for accidental behavior changes, hidden schema changes, ranking changes, and inadequate test coverage.

## Inputs and scope

- Read `.conductor/context-refactor-stage0.md` first.
- Review only the refactor diff relative to baseline `96787c57a0e18baf69f4e48f59be2e7316afffbd` and the relevant tests.
- Do not edit production or test code.

## Checks

- Reference-first five-tier priority remains unchanged.
- QC remains audit-only unless explicit strict override is enabled.
- Prokaryote and eukaryote paths remain distinct where intended.
- Audit, shortlist, review, and selected row schemas and ranking fields remain unchanged.
- Exact accession matching and BUSCO validation invariants remain unchanged.
- Empty, one-candidate, paired GCA/GCF, and override-block paths remain covered.

## Deliverable

Write `.conductor/handoffs/T006-handoff.md` using the conductor handoff contract. Report Critical/Important/Minor findings and an approve/request-changes verdict.

## Claim boundary

You may assess the diff and tests. You must not claim final verification because the conductor owns the final machine-check rerun.
