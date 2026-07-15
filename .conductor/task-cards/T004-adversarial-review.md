# Task card: T004 - Adversarial release review

## Goal
Try to falsify completion of the CDS annotation/strand consistency QC before the
conductor promotes the results. Review implementation, tests, generated summaries,
input manifests, report claims, and Git state.

## Authority / inputs
- `.conductor/context-packet.md`
- `../../../TASK_CDS_ORTHOLOGY_STRAND_CHECK.md`
- `src/align/maf_io.py`, `src/analysis/*.py`, `test/test_*maf*`, `test/test_cds_*`
- `config/cds_orth_strand_inputs.tsv`, `config/ncbi_annotation_availability.tsv`
- `doc/PRELIMINARY_CDS_ORTHOLOGY_CHECK.md`, `doc/results/*.csv`
- External package `data/results/tables/cds_orth_strand_*.csv`

## CLAIM BOUNDARY (hard)
- Review only. Do not modify production, tests, config, docs, results, or Git history.
- You may write only `.conductor/handoffs/T004-handoff.md`.
- Treat every completion claim as unproven until checked against artifacts.
- Flag any scientific overclaim, silent drop, schema mismatch, stale number,
  performance/memory issue, or missing test as a release blocker or residual risk.

## Steps
1. Compare task requirements and corrected context with implementation and outputs.
2. Inspect high-risk coordinate, overlap, classification, and streaming paths.
3. Recount tracked summary/sensitivity values against external generated CSVs.
4. Run focused tests and `git diff --check`.
5. Write a compact handoff with findings ordered by severity.

## Machine checks
- `python3 -m unittest test.test_maf_io test.test_maf_cut_region test.test_cds_intervals test.test_cds_orth_strand_check`
- `cmp doc/results/cds_orth_strand_summary.csv ../data/results/tables/cds_orth_strand_summary.csv`
- `cmp doc/results/cds_orth_strand_sensitivity.csv ../data/results/tables/cds_orth_strand_sensitivity.csv`
- `git diff --check`

## Stop conditions
- Do not repair findings. Report them with paths and exact impact.
