# Handoff: T004 - Adversarial release review blocked promotion

task_id: T004
status: blocked
claim_boundary_ok: true

## Findings

### BLOCKER 1 - Versioned result copies are ignored and absent from Git

`doc/results/cds_orth_strand_summary.csv` and
`doc/results/cds_orth_strand_sensitivity.csv` exist locally and compare
byte-for-byte equal to the external package outputs, but neither file is
tracked or visible as an untracked candidate. `.gitignore:2` ignores every
`results/` directory. A commit made from the current status would therefore
omit the two result tables that T004 describes as generated tracked summaries.
The report would refer to counts without preserving the corresponding compact
CSV evidence on the branch. Promotion must wait until the conductor either
adds an explicit ignore exception or intentionally force-adds and verifies
both files.

Evidence:

- `git ls-files --stage -- doc/results` returned no entries.
- `git status --short --ignored -- doc/results` returned `!! doc/results/`.
- `git check-ignore -v` identified `.gitignore:2:results/` for both CSVs.

### BLOCKER 2 - The binding full-discovery machine check is red

`python3 -m unittest discover -s test -p 'test_*.py'` ran 69 tests and failed
with four errors. `test_chi2_context` cannot import `chi2_context`.
`test_single_sbst_2TSVs`, `test_single_sbst_2TSVs_errprob`, and
`test_subRatio` cannot open fixture paths relative to the repository root.
These appear to be legacy test-infrastructure failures rather than regressions
in the reviewed CDS code, but the context packet explicitly requires this
command to be green and makes a machine-check failure a hard stop. Release
therefore needs either repairs followed by a green rerun or an explicit,
recorded baseline waiver by the conductor.

## Checks completed before the hard stop

- Focused suite: 55 tests passed.
- `python3 -m py_compile` for both changed alignment modules and both analysis
  modules: PASS.
- `git diff --check`: PASS.
- Both local `doc/results` CSVs match their external package counterparts with
  `cmp`: PASS.
- The repository and package copies of
  `PRELIMINARY_CDS_ORTHOLOGY_CHECK.md` match byte-for-byte; SHA-256 is
  `cb23d94a206c96f8d161a15b9fb64af879b96dfd241e10f030c5299cef04deb2`.
- Code inspection found the documented minus-strand coordinate formula,
  disjoint CDS strand-state segments, explicit three-row enforcement,
  per-row streaming output, and atomic replacement of each CSV.

## Non-blocking caveat

The T004 card points to `../../TASK_CDS_ORTHOLOGY_STRAND_CHECK.md`, which does
not exist from the repository root. The actual task file is three levels up at
`../../../TASK_CDS_ORTHOLOGY_STRAND_CHECK.md`. This does not alter computed
results, but it breaks the stated reviewer authority path and should be fixed
before the review record is finalized.

## Unfinished review surface

Per the context packet hard stop, no additional commands were run after the
full-suite failure. Independent raw-row recounting, sensitivity reclassification,
and a second coordinate/strand spot-check remain uncompleted by T004. Existing
report claims for those checks must not be attributed to this reviewer.

## Next action

Resolve or explicitly waive the full-suite baseline, make the compact result
CSVs genuinely versioned, correct the task-card path, then rerun T004 from the
raw-output recount step before promotion.

delegation_effect: difficulty=hard; agents=0; useful=found two concrete release
gates before promotion; misses=hard stop prevented independent raw recount;
next_mix=one conductor repair pass followed by a fresh independent reviewer
