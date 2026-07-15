# Handoff: T002 - CDS annotation segmentation implemented

task_id: T002
status: partial
claim_boundary_ok: true

## Summary

Implemented the T003-facing CDS interval interface in
`src/analysis/cds_intervals.py`: immutable `CDSSegment(start0, end0, strand)`,
`CDSIndex.from_gff(path, feature_type="CDS")`, `query(seqid, coord)`,
`segments_for(seqid)`, and the deterministic `seqids` property. Plain-text and
gzip GFF3 inputs are supported. Selected features are converted from one-based
inclusive coordinates to zero-based half-open coordinates, while invalid
coordinates, missing seqids, and strand values other than `+` or `-` fail with
path and line context.

The index uses per-seqid endpoint events to construct sorted, disjoint genomic
segments. Same-strand overlaps, duplicate records, and touching same-strand
features collapse into a single union segment. Opposite-strand overlap becomes
strand `.`, while opposite-strand features that only touch remain separate.
Point lookup uses binary search over precomputed segment starts; absent seqids
return an empty segment tuple and no query match.

TDD order was preserved. The complete focused test file was added first, and
the focused suite was observed failing because the production module did not
exist. Production code was then added and the focused suite passed all 11
tests. The required full discovery command remains red because four pre-existing
tests outside T002 fail: `test_chi2_context` cannot import `chi2_context`, and
three legacy tests cannot find their relative fixture paths. No files outside
the T002-owned implementation, test, and required handoff were modified.

## Changed files

- `src/analysis/cds_intervals.py`
- `test/test_cds_intervals.py`

## Machine checks

- TDD focused RED observation: PASS
- `python3 -m unittest test.test_cds_intervals`: PASS
- `python3 -m unittest discover -s test -p 'test_*.py'`: FAIL
- `python3 -m py_compile src/analysis/cds_intervals.py`: PASS
- `git diff --check -- src/analysis/cds_intervals.py test/test_cds_intervals.py`: PASS

## Unresolved risk / escalation

- Full discovery has four errors in non-owned legacy tests: one import-path
  error and three missing relative-fixture errors. T002 tests themselves pass.
- ESCALATE: JUDGMENT_REQUIRED
  question: Should the conductor accept the known baseline failures for T002,
  or assign their repair to the owners of the legacy test infrastructure?
  evidence: `test/test_chi2_context.py`, `test/test_single_sbst_2TSVs.py`,
  `test/test_single_sbst_2TSVs_errprob.py`, and `test/test_subRatio.py`

## Next action

Have the conductor independently inspect the T002 diff, resolve or waive the
non-owned baseline failures, and rerun the required full discovery command.

delegation_effect: difficulty=simple; agents=0; useful=bounded T002 implementation
and test evidence; misses=legacy full-suite baseline remained red; next_mix=one
conductor reviewer plus the legacy-test owner if full-suite repair is required
