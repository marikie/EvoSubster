# Task card: T003 - Streaming classifier, real inputs, and evidence

## Goal
Build the four-class streaming analysis, acquire and validate exact-assembly GFFs,
run all three trios, and generate auditable outputs and documentation.

## Authority / inputs
- Shared context: `.conductor/context-packet.md`
- Main-session ownership: `src/analysis/cds_orth_strand_check.py`,
  `test/test_cds_orth_strand_check.py`, configuration, docs, external input/result files.

## CLAIM BOUNDARY (hard)
- Allowed: classify only coding/strand consistency with explicit insufficient evidence.
- MUST NOT label results as proof of orthology or modify manuscript claims/visualizations.
- All input provenance, skipped-block accounting, and invariants must be recorded.
- Generated outputs must come from the script, never manual CSV editing.
- No commit until full verification and release decision.

## Steps
1. Write failing tests for all four classes, thresholds, ambiguity, strict blocks,
   streaming output, sensitivity summaries, and atomic failure behavior.
2. Implement classifier and CLI using T001/T002 interfaces.
3. Download six exact-assembly GFFs and record URL, bytes, SHA-256, and gzip/GFF checks.
4. Run `lenEdo`, `puffer`, and `styPis`; verify all accounting invariants.
5. Hand-check evidence rows and generate report/manifest.

## Machine checks
- `python3 -m unittest test.test_cds_orth_strand_check`
- `python3 -m unittest discover -s test -p 'test_*.py'`
- Real-data preflight and full run for all three trios.

## Stop conditions
- Stop before publishing results if an input does not match its MAF assembly,
  an invariant fails, or performance/memory limits prevent a complete run.
