# Handoff: T001 - Eupercaria Stage 0 candidate inventory completed

task_id: T001
status: done
claim_boundary_ok: true
o_level: O1 (PROPOSED; inventory facts only)
paired_non_claim: This inventory does not establish that BUSCO changes, or does not change, any Stage 0 winner.

## Summary

The supplied tree `/Users/marikonakagawa/biohazard/TheTreeOfLife/Eupercaria_tree.nwk` contains 193 tips representing 148 species under the exact `species_from_label()` parsing rule in `src/select/trio_selection.R`. Thirty-three species occur on multiple tips (45 extra duplicate-species tips; the maximum is seven tips for *Larimichthys crocea*). The retained Stage 0 snapshot is under `results/eupercaria_20260801/trio_selection/`: `assembly_metadata.tsv` has 375 records for 145 species, `assembly_candidates_audit.tsv` has the same 375 records (365 eligible, 10 excluded), `assembly_shortlist.tsv` has 305 rows for 145 species, and `selected_assemblies.tsv` has one winner for each of those 145 species. The three parsed names absent from metadata are `Chaunax sp.`, `Pseudoliparis sp.`, and `Platycephalus sp.`.

All 145 baseline winners are marked `is_reference=TRUE`; 47 are annotated GCF winners and 77 winners have some annotation. One hundred species have at least two eligible/shortlisted candidates, so real alternatives exist. Baseline selection uses `reference_metadata_provisional` for 52 species and `reference_metadata_general` for 93. External `qc_busco_*` fields are empty throughout this snapshot, so it is a valid Reference-first baseline but contains no evidence about a BUSCO override.

The most bounded scientifically useful experiment is *Takifugu rubripes*, which is present in the supplied tree and has an existing focused three-species result at `results/eupercaria_20260801/takifugu_e2e/`. Its baseline winner is the annotated RefSeq Reference `GCF_901000725.3`. Its current eligible shortlist also contains `GCA_901000725.3` (the GenBank counterpart with the same accession stem and identical recorded assembly statistics) and the distinct non-reference haploid Chromosome assembly `GCA_000180615.2`. The reference FASTA already exists at `genomes/eupercaria_20260801/Takifugu_rubripes/GCF_901000725.3_fTakRub1.3_genomic.fna`; no FASTA for `GCA_000180615.2` was found. For a fair bounded BUSCO test, run the same genome-mode BUSCO version, exact Actinopterygii lineage dataset/version, and settings on the reference and `GCA_000180615.2`. Avoid treating `GCA_901000725.3` as an independent biological assembly; verify FASTA checksums before reusing the GCF BUSCO score for that paired record.

`GCA_901000745.3` is a current `alternate-pseudohaplotype` and has shortlist rank 4, so it is outside the default top three. The current unsupported-type regex does not exclude this spelling, although project prose says alternate haplotypes are hard-filtered. It should not enter the BUSCO comparison without a conductor decision.

A command that regenerates the Stage 0/current-NCBI snapshot and dry-run Stage 1 layout is:

`Rscript src/select/trio_selection.R --tree /Users/marikonakagawa/biohazard/TheTreeOfLife/Eupercaria_tree.nwk --out-dir results/eupercaria_20260801/trio_selection --date 20260801 --stage0-top-k 3 --dry-run`

Because NCBI metadata is live, this reproduces the procedure, not necessarily the exact August 1 records.

## Changed files

- `.conductor/handoffs/T001-handoff.md`

## Machine checks

- `git status --short --branch` PASS (target branch; only pre-existing/untracked `.conductor/`, `analysis/`, and HTML artifacts)
- supplied Newick exists PASS
- four retained Stage 0 artifact paths exist PASS
- tree/tables parsed and counted PASS

## Unresolved risk / escalation

- ESCALATE: JUDGMENT_REQUIRED
  question: Should the bounded experiment treat the GCA/GCF counterpart as one assembly after checksum verification and compare only `GCF_901000725.3` versus distinct `GCA_000180615.2`, while excluding the alternate pseudohaplotype?
  evidence: `results/eupercaria_20260801/trio_selection/assembly_candidates_audit.tsv`, `src/select/assembly_selection.R`

## Next action

The conductor should verify the Takifugu accession scope, download only the missing distinct assembly, run like-for-like BUSCO, build the exact-accession QC TSV, and independently compare baseline versus QC-ranked selections.

delegation_effect: difficulty=medium; agents=1; useful=identified the smallest non-pseudoreplicated real comparison and the alternate-pseudohaplotype filter gap; misses=no BUSCO execution by design; next_mix=one execution worker plus one independent evidence reviewer
