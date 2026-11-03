# Task card: T003 - Execute the predeclared BUSCO pilot

## Goal
Download the six exact candidate genome FASTAs in the predeclared protocol, run BUSCO genome mode under identical conditions, and produce raw summary artifacts plus a schema-compatible external-QC TSV draft.

## Authority / inputs
- Context: `.conductor/context-packet.md`
- Protocol: `results/eupercaria_20260807/busco_experiment/experiment_protocol.md`
- Candidate list: `results/eupercaria_20260807/busco_experiment/candidate_accessions.tsv`
- Runtime: `results/eupercaria_20260807/busco_experiment/runtime/venv`
- Lineage cache: `results/eupercaria_20260807/busco_experiment/busco_downloads`

## CLAIM BOUNDARY (hard)
- Allowed: download the six listed FASTAs, run the pinned BUSCO recipe, preserve logs/summaries/checksums, and draft the exact-accession TSV.
- MUST NOT: add candidates after seeing results, edit application source, promote an override conclusion, delete user files, commit, or push.
- The paired GCA Reference row may reuse the GCF result only as explicitly declared in the protocol, and provenance must remain visible.

## Steps
1. Read the context packet, protocol, and candidate list.
2. Verify the BUSCO, Miniprot, HMMER, and amended `actinopterygii_odb10` lineage versions/checksums.
3. Download only genome FASTA for all six listed input accessions; record SHA-256 and file sizes.
4. Run every input with BUSCO 6.1.0, `-m genome`, explicit local `actinopterygii_odb10`, `--offline`, `--miniprot`, `--skip_bbtools`, and `-c 3`.
5. Parse `short_summary` JSON files into `external_busco_qc.tsv`, including paired GCA Reference rows with `qc_input_accession` provenance.
6. Report failures honestly and stop rather than substituting a different tool or dataset.

## Machine checks
- All six FASTA files are non-empty and begin with `>`.
- All six BUSCO runs exit successfully and have parseable `short_summary` JSON.
- TSV has exact versioned accessions, mode `genome`, one lineage/version, numeric Complete/Fragmented/Missing values, and no duplicate accessions.
- `git status --short --branch`.

## Stop conditions
- On a dependency or reproducibility failure, stop and report it; do not switch BUSCO versions or lineage datasets.
- If disk free space falls below 10 GiB, stop before launching another run.

## Handoff
Write `.conductor/handoffs/T003-handoff.md` (400-800 words, paths and pass/fail lines only; no raw logs).
