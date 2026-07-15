# Decision: D-000 - Waive one unrelated stale-import baseline failure

Decision: The CDS annotation/strand QC release gate may proceed when its focused
suite is green even though repository-wide discovery retains one pre-existing
`test_chi2_context` import error. This is a narrow baseline waiver, not a claim
that the full suite passes.

Evidence:

- From `test/`, `python3 -m unittest discover -p 'test_*.py'` runs 70 tests;
  only `test_chi2_context` errors because `chi2_context` cannot be imported.
- Commit `b6b33ae60fb4bf350fbcd60e042a039d37fa8269` intentionally deleted
  `src/statistics/chi2_context.py` as superseded, while leaving the stale test.
- The current change does not modify `src/statistics/`, `test/test_chi2_context.py`,
  or the legacy fixture-based tests.
- The four apparent errors from repository-root discovery include three
  working-directory-sensitive legacy fixture paths. Running discovery from
  `test/`, as those tests expect, removes those three errors.
- The task-focused command runs 56 tests and passes.

Conditions:

- T004 must continue its remaining raw-row recount, sensitivity, coordinate,
  strand, report, and Git-state checks after reading this decision.
- Any new failure in a task-owned test remains a blocker.
- The final release decision must disclose the one waived baseline error and
  must not describe repository-wide discovery as green.

Paired non-claim: This waiver does not validate or restore the removed
`chi2_context.py` workflow and does not excuse failures caused by this task.
