# Task card: T004 - Independently verify and try to refute the BUSCO override

## Goal
Re-derive all six BUSCO scores directly from the six specific JSON summaries, compare them against the worker TSV, run the current Stage 0 ranking on the saved metadata, and try to refute any apparent Reference override.

## Authority / inputs
- Context: `.conductor/context-packet.md`
- Protocol and candidates: `results/eupercaria_20260807/busco_experiment/experiment_protocol.md`, `candidate_accessions.tsv`
- Raw summaries: `results/eupercaria_20260807/busco_experiment/busco_runs/*/short_summary.specific.*.json`
- Draft TSV: `results/eupercaria_20260807/busco_experiment/external_busco_qc.tsv`
- Baseline: `results/eupercaria_20260801/trio_selection/assembly_metadata.tsv`
- Ranking code: `src/select/assembly_selection.R`

## CLAIM BOUNDARY (hard)
- Allowed: create an independent derived TSV/comparison under the ignored experiment directory, report exact discrepancies, and propose accept/reject/downgrade for the bounded claim.
- MUST NOT: edit the worker TSV, raw summaries, application source, candidate scope, or protocol; commit or push; claim anything about the other 142 retained species.
- Treat GCF-to-paired-GCA changes as the same Reference assembly, not a non-reference override.
- Try to refute: check score units, JSON field mapping, accession duplication/provenance, full shortlist coverage, ranking tie-breaks, and internal-stop warnings.

## Steps
1. Parse the six JSONs independently and confirm mode/lineage/software/input identity.
2. Reconstruct all nine exact-accession rows and byte-compare the score fields with the draft TSV.
3. Confirm all saved top-three shortlist rows for the three species have external QC.
4. Run baseline and QC Stage 0 rankings using current code; identify true non-reference overrides separately from paired-accession changes.
5. Inspect whether any warning or numerical tie invalidates the claim.
6. Write and lint a compact adversarial handoff.

## Machine checks
- Six unique measured accessions and nine unique TSV accessions.
- Complete + Fragmented + Missing equals 100 within summary rounding tolerance for every measured run.
- One mode, lineage, and BUSCO version across rows.
- Current R Stage 0 tests pass and the saved snapshot rerank completes.
- `git status --short --branch`.

## Stop conditions
- Any raw-summary/TSV mismatch: stop and recommend rejection until corrected.
- Any missing shortlist accession: do not accept a quality-first claim for that species.

## Handoff
Write `.conductor/handoffs/T004-handoff.md` (400-800 words, concise evidence paths and pass/fail checks, proposed O-level and paired non-claim).
