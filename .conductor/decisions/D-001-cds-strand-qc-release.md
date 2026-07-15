# D-001 - Accept bounded three-trio CDS strand QC

Date: 2026-07-15
Tier: 0 (independent Codex judgment panel)
Conductor: Codex main session
Inputs judged: `.conductor/handoffs/T001-handoff.md`, `T002-handoff.md`,
`T003-handoff.md`, `T004-handoff.md`, `T005-handoff.md`, and
`.conductor/decisions/D-000-existing-test-baseline.md`

## Question

Can the implementation and generated evidence be accepted as a reproducible,
bounded quality-control result for independent CDS coverage and relative strand
consistency in the exact `lenEdo`, `puffer`, and `styPis` inputs? An incorrect
acceptance would make later visualization decisions depend on coordinate,
annotation, or classification artifacts.

## Verdict

VALID within the documented claim boundary. The implementation, compact
summaries, external raw evidence, provenance manifest, and report agree. T004's
two initial blockers were repaired, and T005 independently passed the repaired
state. D-000 explicitly waives only one unrelated stale-import baseline error.

## Outcome level + paired non-claim

- Ledger transition: draft -> reviewed -> accepted
- O-level: O3, because this is an evidence-bounded exploratory QC with complete
  accounting and sampled independent coordinate validation, but not external
  biological or orthology validation.
- Paired non-claim: The accepted result does not prove orthology, conserved
  synteny, absence of paralogy, annotation correctness, or a biological cause
  for coding-state differences.
- Independent judgment panel: agreed after the first reviewer blocked promotion;
  the conductor repaired Git tracking and CRLF generation, recorded the narrow
  test-baseline waiver, and a fresh reviewer passed the final state.

## Load-bearing reasoning

All 449,420 raw rows were independently recounted, and their primary class
distributions match the generated summary. Independent reclassification from
raw evidence reproduced sensitivity thresholds 1, 10, and 30. Nine rows spanning
three substantive classes and all three trios were rebuilt directly from MAF
and GFF inputs, including minus-strand forward coordinates; reported coverage,
strand, and class values matched. The six GFF inputs have versioned exact
assembly accessions, byte counts, SHA-256 values, CDS presence, and MAF seqid
compatibility checks. Large inputs and raw evidence remain external; the branch
tracks the implementation, provenance, report, and compact result tables.

## Conditions attached

- Preserve the exact input checksums and corrected puffer B accession
  `GCF_003711565.1` when reproducing these counts.
- Keep `insufficientEvidence` separate and retain the support threshold in any
  downstream visualization.
- Do not rename `codingStrandConsistent` to an orthology claim.
- Report repository-wide tests honestly: 56 task tests pass; full discovery
  runs 70 tests with the one D-000 stale `chi2_context` import error.

## Applied

The streaming parser, CDS interval index, classifier, CLI, tests, provenance
tables, report, and compact summaries were generated or reviewed in the paths
listed by T003. Primary counts are:

| trio | non-coding | consistent | contradictions | insufficient | total |
|---|---:|---:|---:|---:|---:|
| lenEdo | 10,283 | 57,972 | 98 | 2,592 | 70,945 |
| puffer | 42,368 | 175,159 | 25 | 4,091 | 221,643 |
| styPis | 15,547 | 135,927 | 44 | 5,314 | 156,832 |

Focused tests, byte compilation, manifest preflight, full real-data execution,
raw accounting, sampled direct reconstruction, report/CSV comparisons, LF
regression, and staged whitespace checks pass. The combined run completed in
812.1 seconds at about 270 MB peak RSS.

## Residual / follow-up

The stale `test_chi2_context` baseline and non-fatal `test_subRatio` resource
warning remain outside this task. Any manuscript-facing interpretation requires
separate orthology/synteny analysis and domain review. Exhaustive raw accounting
was independent, but direct MAF/GFF coordinate reconstruction was sampled rather
than exhaustive.

delegation_effect: difficulty=hard; agents=4 across implementation and two
review passes; useful=bounded modules, compatibility evidence, two concrete
release-blocker discoveries, and independent raw/coordinate verification;
misses=the first review stopped before raw recount and CRLF appeared only after
the ignored CSVs were staged; next_mix=use two bounded implementation workers,
stage compact artifacts before review, then run one fresh read-only gate.
