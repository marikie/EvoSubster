# Handoff: T005 - Independent release re-review passed promotion

Status: done
claim_boundary_ok: true
o_level: O3 (PROPOSED; the conductor decides; this review does not self-promote the claim)
paired_non_claim: This QC does not prove orthology, conserved synteny, paralogy status, annotation correctness, or a biological cause for coding-state differences.

## Summary

The final workspace state is suitable for promotion within the stated CDS
annotation and relative-strand consistency claim boundary. No release blocker
remains. The two T004 gates were addressed: both compact result CSVs are now
present in the Git index, and repository-wide discovery from `test/` differs
from green only by the exact stale-import baseline waived in
`.conductor/decisions/D-000-existing-test-baseline.md`.

During this review, an intermediate staged state exposed CRLF line endings in
the compact CSVs through `git diff --cached --check`. The conductor corrected
the production CSV writer, added a regression test, regenerated both external
and repository copies with LF, and re-staged them. The final review reran all
affected gates: focused discovery now passes 56 tests, both unstaged and staged
whitespace checks pass, and each staged compact CSV is byte-identical to its
external generated counterpart. The report copies also remain byte-identical.

An independent raw recount read all 449,420 rows without importing the analysis
module. Raw class totals match each threshold-10 summary row: `lenEdo` has
10,283 non-coding, 57,972 consistent, 98 contradictions, and 2,592
insufficient; `puffer` has 42,368, 175,159, 25, and 4,091; `styPis` has 15,547,
135,927, 44, and 5,314. A separate reclassification from raw evidence reproduced
every class count and accounting invariant at thresholds 1, 10, and 30.

Nine representative rows were then reconstructed directly from the compressed
MAFs and six GFFs: one `codingStrandConsistent`, one `nonCodingIncluded`, and one
`strandContradictions` row per trio. Forward coordinates were recalculated for
all three MAF rows, including minus rows with the reverse-coordinate formula.
Direct GFF interval scans reproduced non-gap counts, CDS coverage, opposite-
strand conflicts, majority strand calls, and final classes for all nine rows.

Report numbers, the corrected `../../../TASK_CDS_ORTHOLOGY_STRAND_CHECK.md`
authority paths in T004/T005, and required artifact paths were checked. No file
larger than 5 MB is present outside `.git`; large MAF, GFF, and raw result files
remain in the external package as intended.

## Changed files

- `.conductor/handoffs/T005-handoff.md`

No production code, tests, configs, docs, results, Git index, or Git history
were modified by this reviewer.

## Machine checks

- Focused task suite, 56 tests  PASS
- Full discovery from `test/`, 70 tests  FAIL only at waived `test_chi2_context` stale import
- Git index presence for both `doc/results/*.csv` files  PASS
- Report and both compact CSV external/repository byte comparisons  PASS
- Independent 449,420-row recount and threshold 1/10/30 reclassification  PASS
- Nine independent MAF-coordinate and GFF-strand checks  PASS
- `git diff --check` and `git diff --cached --check` on final state  PASS
- Task-card paths and unexpected large repository files  PASS

## Unresolved risk / escalation

- No blocker.
- Repository-wide discovery is not green; D-000 waives only the pre-existing
  removed `chi2_context.py` import. The unrelated `test_subRatio` resource
  warning also remains non-fatal.
- Reproduction depends on external large inputs and raw outputs identified by
  the versioned manifest and checksums.
- Coordinate validation is a nine-row independent sample; exhaustive raw
  accounting was performed, but exhaustive independent GFF reconstruction was
  not.

## Next action

The conductor may record the release decision and commit the reviewed artifacts
on `feature/protein-coding-region`, while retaining D-000 and the paired
non-claim in the release record.

delegation_effect: difficulty=hard; agents=0; useful=independent recount and
direct MAF/GFF reconstruction closed the T004 unfinished review surface;
misses=an intermediate CRLF defect required conductor repair and a final rerun;
next_mix=retain a read-only reviewer after integration and staging are complete.
