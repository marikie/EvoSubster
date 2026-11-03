# D-001 - Accept current-run unused species cleanup

Date: 2026-08-08
Tier: other - hard destructive-default code gate with independent adversarial review
Conductor: Codex (GPT-5)
Inputs judged: `.conductor/handoffs/T001-handoff.md`, `.conductor/handoffs/T002-handoff.md`, `.conductor/handoffs/T003-handoff.md`, `docs/superpowers/specs/2026-08-08-unused-trio-species-cleanup-design.md`

## Question

Does the implementation safely delete only current-run genome and LAST
artifacts for accessions absent from every selected trio, while preserving
selected and pre-existing data and stopping downstream analysis on cleanup
failure?

## Verdict

VALID within the fixture-tested contract. The change is accepted for a local
commit after the downloader manifest, selector cleanup, wrapper propagation,
zero-trio behavior, audit failure reporting, and LAST database ownership all
survived main-session reproduction and an independent adversarial pass.

## Outcome level + paired non-claim

- Ledger transition: draft -> reviewed -> accepted
- O-level: O3 - evidence-bounded implementation with complete local fixture and regression coverage, but no live NCBI/LAST run.
- Paired non-claim: Passing local tests and adversarial probes does not establish live NCBI/LAST behavior or safety under concurrent processes mutating the same genome directory.
- Independent judgment panel: T003 found four concrete gaps across two review rounds; all were converted to regressions, fixed, and rechecked successfully.

## Load-bearing reasoning

The downloader reports only paths absent before its successful call and keeps
its six-field stdout stable. The selector accepts only root-relative manifest
paths, rejects traversal and parent-symlink escapes before deletion, does not
follow direct symlinks, and treats recursive ownership conflicts as fatal. It
retains genome artifacts by selected accession, train files by selected pair,
and LAST databases by selected reference accession. The audit is atomically
replaced before preflight and records the active failing row. Wrapper tests show
that selector failure launches no downstream run and that a blank zero-trio
table exits successfully.

## Conditions attached

- Keep `--keep-unused-species-data` as the explicit debugging opt-out.
- Do not claim concurrent shared-directory safety.
- A future live smoke test must use disposable output directories and must not
  broaden the deletion contract without another destructive-path review.

## Applied

Changed the downloader, selector, wrapper, README, design document, and focused
tests. Main verification passed all three shell integration suites, all three R
suites, 44 Python tests, syntax checks for every shell script, and
`git diff --check`. The Python suite still emits its pre-existing unclosed-file
`ResourceWarning`; it does not fail.

## Residual / follow-up

Live NCBI and LAST execution was not performed because the accepted test scope
uses deterministic process-boundary fixtures. No release blocker remains in
that scope.

delegation_effect: difficulty=hard; agents=3; useful=parallel implementation isolated write scopes and the adversarial reviewer found four actionable gaps; misses=live external-tool behavior; next_mix=two implementers plus one focused adversarial reviewer for similar destructive pipeline changes.
