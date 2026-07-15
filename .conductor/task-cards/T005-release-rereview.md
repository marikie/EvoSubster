# Task card: T005 - Release re-review after gate repairs

## Goal
Independently determine whether the CDS annotation/strand consistency QC can be
promoted after the first adversarial review's actionable gates were addressed.

## Authority / inputs
- `.conductor/context-packet.md`
- `.conductor/handoffs/T004-handoff.md`
- `.conductor/decisions/D-000-existing-test-baseline.md`
- `../../../TASK_CDS_ORTHOLOGY_STRAND_CHECK.md`
- Implementation, tests, configs, tracked report/summary CSVs, and external raw
  outputs named in T004.

## CLAIM BOUNDARY (hard)
- Review only; write only `.conductor/handoffs/T005-handoff.md`.
- The D-000 waiver is limited to the one documented pre-existing stale import.
- Do not attribute prior independent spot checks to this review. Recheck enough
  raw evidence directly to support or reject promotion.
- This QC does not prove orthology, synteny, paralogy status, or biological cause.

## Steps
1. Verify both `doc/results` CSVs are now present in the Git index and still
   byte-identical to external generated outputs.
2. Run the focused task suite and the full suite from `test/`; confirm the
   latter differs from green only by the exact D-000 baseline.
3. Independently recount raw row/class totals and summary/sensitivity invariants.
4. Independently inspect representative minus-coordinate and strand-class rows.
5. Check report claims, task-card paths, whitespace, and unexpected large Git files.
6. Issue a promotion/pass or blocker finding with residual risks.

## Machine checks
- `python3 -m unittest test.test_maf_io test.test_maf_cut_region test.test_cds_intervals test.test_cds_orth_strand_check`
- From `test/`: `python3 -m unittest discover -p 'test_*.py'`
- `git ls-files --stage -- doc/results/*.csv`
- `cmp` for report and both summary CSV copies
- `git diff --check`

## Stop conditions
- Do not repair findings. Report them with paths and exact impact.
