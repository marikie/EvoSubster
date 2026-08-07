# T003 handoff

task_id: T003-run-busco-pilot

status: done

## summary

The six predeclared NCBI genome packages were downloaded without adding or substituting candidates. Package files and transfer logs are under `results/eupercaria_20260807/busco_experiment/downloads/`. All six ZIP archives are structurally valid. The extracted genome FASTAs are `results/eupercaria_20260807/busco_experiment/inputs/<ACCESSION>.fna`; every file is non-empty and begins with `>`. Exact byte counts and SHA-256 values are recorded in `results/eupercaria_20260807/busco_experiment/input_manifest.tsv`.

The lineage archive `results/eupercaria_20260807/busco_experiment/busco_downloads/lineages/actinopterygii_odb10.2024-01-08.tar.gz` was verified at 177,093,579 bytes with MD5 `6b042870822559fcacf129fdedea8927`. Its extracted dataset is `results/eupercaria_20260807/busco_experiment/busco_downloads/lineages/actinopterygii_odb10`; `dataset.cfg` records creation date `2024-01-08` and 3,640 BUSCO groups. Runtime evidence records BUSCO 6.1.0, Miniprot 0.18-r281, HMMER 3.4, and the conductor-provisioned BBTools 40.02 dependency.

An initial configuration attempt for `GCF_901000725.3` failed because BUSCO checked for BBTools despite `--skip_bbtools`. That diagnostic remains preserved at `results/eupercaria_20260807/busco_experiment/busco_logs/GCF_901000725.3.initial_failure.log`. After the conductor installed BBTools 40.02, the failed output was replaced with `--force` and the same pinned biological conditions. No candidate, BUSCO version, lineage, predictor, scoring rule, or CPU setting changed.

All six final runs completed with genome mode, the explicit local `actinopterygii_odb10` lineage, offline mode, Miniprot, `--skip_bbtools`, and three CPUs. Each run has one specific parseable JSON at `results/eupercaria_20260807/busco_experiment/busco_runs/<ACCESSION>/short_summary.specific.actinopterygii_odb10.<ACCESSION>.json`; complete, fragmented, and missing percentages are numeric. Final run logs are under `results/eupercaria_20260807/busco_experiment/busco_logs/`. Every run emitted a Miniprot internal-stop-codon warning, ranging from 6.0% to 9.2% of Complete matches; these warnings are preserved and are not silently recoded.

The schema-compatible draft is `results/eupercaria_20260807/busco_experiment/external_busco_qc.tsv`. It contains nine unique exact versioned accessions: six measured input rows plus three paired GCA Reference rows. Paired rows reuse their declared GCF measurement and preserve the measured accession in `qc_input_accession`. All rows record mode `genome`, lineage `actinopterygii_odb10`, BUSCO version `6.1.0`, and numeric Complete/Fragmented/Missing values. This handoff reports execution evidence only and does not promote or interpret a Stage 0 Reference-override claim.

## changed_files

- `.conductor/handoffs/T003-handoff.md`
- `results/eupercaria_20260807/busco_experiment/input_manifest.tsv`
- `results/eupercaria_20260807/busco_experiment/inputs/`
- `results/eupercaria_20260807/busco_experiment/downloads/`
- `results/eupercaria_20260807/busco_experiment/busco_logs/`
- `results/eupercaria_20260807/busco_experiment/busco_runs/`
- `results/eupercaria_20260807/busco_experiment/external_busco_qc.tsv`

## machine_checks

- PASS: all six exact package ZIPs are valid; all six FASTAs are non-empty and FASTA-prefixed.
- PASS: the input manifest contains six exact versioned accessions with byte counts and SHA-256 values.
- PASS: odb10 archive size, published MD5, dataset date, and 3,640-group count match the amended protocol.
- PASS: all six final logs report BUSCO completion and contain no `ERROR:` line.
- PASS: exactly six specific summary JSON files exist and validate identical pinned parameters and software versions.
- PASS: `external_busco_qc.tsv` has nine data rows, no duplicate accession, one mode, one lineage, one BUSCO version, and numeric score fields.
- PASS: the TSV passed current Stage 0 attachment and comparable-BUSCO validation.
- PASS: no BUSCO, Miniprot, or HMMER process remained after completion; free disk stayed above the 10 GiB stop threshold.
- PASS: `git status --short --branch` confirmed the requested feature branch; no commit or push was performed.

claim_boundary_ok: true

next_action: The conductor should independently verify the six summaries and run the TSV through the saved Stage 0 shortlist before accepting, downgrading, or rejecting any Reference-override claim.

delegation_effect: difficulty=hard; agents=1; useful=produced six checksummed inputs, six same-condition BUSCO summaries, a provenance-preserving nine-row TSV, and exposed then recovered the hidden BBTools startup dependency before results were observed; misses=no release-gate interpretation was attempted by design; next_mix=one independent result/TSV verifier plus conductor judgment is sufficient.
