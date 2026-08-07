# T004 handoff

task_id: T004-adversarial-verify

status: done

## summary

I independently re-derived all six BUSCO measurements from the specific JSON summaries, without using the draft QC TSV as the source of scores. The derived nine-accession table is `results/eupercaria_20260807/busco_experiment/adversarial_busco_qc.tsv`; three paired GCA Reference rows reuse the declared GCF measurement while retaining `qc_input_accession`. A field-level comparison of accession, provenance, mode, lineage, version, Complete, Fragmented, and Missing against `external_busco_qc.tsv` found no mismatch. All nine saved top-three shortlist accessions for the three pilot species are covered.

The six JSONs agree on BUSCO 6.1.0, BUSCO internal mode `euk_genome_min` (the genome-mode implementation), Miniprot 0.18-r281, HMMER 3.4, three CPUs, `actinopterygii_odb10` dated 2024-01-08, and 3,640 markers. Each JSON input path ends in its exact measured accession FASTA. For every run, Complete + Fragmented + Missing is 100.0% after one-decimal reporting, and the corresponding integer counts sum exactly to 3,640. NCBI package `assembly_data_report.jsonl` independently confirms the three GCF/GCA pairings used for score propagation.

Running the current `rank_assembly_candidates()` on the saved metadata and the independently derived QC table reproduced the worker result. `Takifugu rubripes` changes from `GCF_901000725.3` to `GCA_901000725.3`, but these are the paired accessions of the same Reference assembly and have identical BUSCO scores; the change is an accession tie-break, not a non-reference override. `Plectropomus leopardus` changes from Reference `GCF_008729295.1` to distinct non-reference `GCA_026936395.1`: Complete rises from 3,502/3,640 (96.2%) to 3,635/3,640 (99.9%), Fragmented falls from 77 to 0, Missing from 61 to 5, and Duplicated from 46 to 11. This survives adversarial review as a clear bounded BUSCO-driven override under the current Stage 0 rule.

`Larimichthys crocea` mechanically changes from Reference `GCF_000972845.2` to distinct non-reference `GCA_059630365.1`, because Complete rises by four markers (99.7% to 99.8%). I recommend not promoting this as a clearly better assembly: Missing worsens from 3 to 6 and Duplicated rises from 10 (0.3%) to 215 (5.9%), a dimension omitted from the external TSV/ranking. The exact adversarial comparison is recorded in `results/eupercaria_20260807/busco_experiment/adversarial_rerank_comparison.tsv`.

All final runs warn that 6.0-9.2% of Complete Miniprot predictions contain internal stop codons. The warning affects every candidate under the same predictor and does not erase the Plectropomus separation, but it limits interpretation beyond this ranking experiment. The current R metadata also represents blank non-reference categories as `is_reference=NA`, not `FALSE`; exact NCBI accessions and categories establish the two distinct non-reference winners, and this NA does not change these winners because BUSCO fields separate them earlier in the ordering.

## changed_files

- `.conductor/handoffs/T004-handoff.md`
- `results/eupercaria_20260807/busco_experiment/adversarial_busco_qc.tsv`
- `results/eupercaria_20260807/busco_experiment/adversarial_rerank_comparison.tsv`

## machine_checks

- PASS: six unique measured accessions; nine unique QC accessions; six unique provenance inputs.
- PASS: raw JSON fields match every independently reconstructed row and the draft TSV score fields.
- PASS: all nine saved top-three shortlist rows have external QC.
- PASS: one BUSCO version, internal genome mode, lineage dataset/version, predictor, and materially identical settings across six runs.
- PASS: percentage totals are 100.0 within reporting precision; integer totals are exactly 3,640.
- PASS: independent baseline/QC snapshot rerank completed and distinguished paired-accession changes from true overrides.
- PASS: `test_stage0_gate.R`, `test_stage0_quality_selection.R`, and `test_trio_selection_helpers.R` passed.
- PASS: required artifacts exist and `git status --short --branch` remains on the requested feature branch; no commit or push.

claim_boundary_ok: true

o_level: O3 (proposed) - Accept the bounded claim that current Stage 0 can replace a Reference with a distinct non-reference in this predeclared three-species Eupercaria pilot; Plectropomus is the clear positive case.

paired_non_claim: This does not establish that the non-reference is globally superior for every downstream use, that Larimichthys is clearly improved, or that any of the other 142 retained Eupercaria species would be overridden.

next_action: The conductor should promote only the bounded Plectropomus override, describe Larimichthys as a mechanical but duplication-confounded override, and explicitly reject the Takifugu paired-accession change as evidence of a Reference override.

delegation_effect: difficulty=hard; agents=1; useful=separated one false paired-accession change from two algorithmic overrides and downgraded the duplication-confounded Larimichthys interpretation while preserving the clear Plectropomus case; misses=no independent rerun with another BUSCO lineage or gene predictor was in scope; next_mix=the current conductor plus one adversarial verifier is sufficient for this bounded release gate.
