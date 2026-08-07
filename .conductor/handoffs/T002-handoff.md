# T002 handoff

task_id: T002-busco-runtime

status: done

## summary

The current Mac cannot run a reproducible BUSCO comparison without first adding an isolated runtime. It is macOS 26.5.2 on arm64 with 10 logical CPUs, 24 GiB RAM, and 51 GiB free on the data volume. `busco`, Conda/Mamba/Micromamba, Docker/Podman, Miniprot, MetaEuk, Augustus, Prodigal, HMMER, BLAST, and BBMap are absent from `PATH`; `/usr/bin/python3` is 3.9.6 and has no `busco` module. Homebrew 6.0.12 and Java 24.0.2 are present, but neither supplies a complete BUSCO pipeline. No local BUSCO download cache or common Conda installation directory exists. Because upstream BUSCO explicitly does not support macOS and recommends its container on incompatible systems, the least risky reproducible route is the official BUSCO Linux container on a Docker-capable Linux host, not a manual native-Mac installation.

The appropriate common dataset for the Eupercaria fish assemblies is `actinopterygii_odb12.2`, explicitly selected rather than auto-lineage. The live official `file_versions.tsv` reports creation date `2026-05-13` and MD5 `4b3f2a3481ae6813f04e0896dad9b805`; the current archive is 1,481,330,473 bytes compressed. BUSCO's current guide says odb12.2 requires BUSCO 6.1.0 or later, identifies Miniprot as the default eukaryotic genome predictor, recommends the most specific fish lineage (`actinopterygii`), and requires reporting dataset creation date. Sources: <https://busco.ezlab.org/busco_userguide.html> and <https://busco-data.ezlab.org/v5/data/file_versions.tsv>.

Recommended pinned recipe: use `ezlabgva/busco:v6.1.0_cv2`; download `actinopterygii_odb12.2` once into a mounted dataset directory; verify the recorded archive MD5/date; then run every nucleotide FASTA with `-m genome`, the same absolute local lineage path, `--offline`, explicit `--miniprot`, the same `-c` value, default scoring parameters, `--tar`, and separate output names. Record container tag/digest, command, input FASTA checksum, exact versioned accession, dataset name/date, and the `short_summary` JSON. Convert its percentages into the repository schema columns `accession`, `qc_busco_mode=genome`, `qc_busco_lineage=actinopterygii_odb12.2`, `qc_busco_version=6.1.0`, `qc_busco_complete`, `qc_busco_fragmented`, and `qc_busco_missing`. Do not mix NCBI annotation-protein BUSCO values with these genome-mode results.

The saved Eupercaria shortlist has 305 accession rows across 145 species, representing 258 unique GCA/GCF assembly keys and 194.53 Gbp of assembled sequence; 92 species have both a Reference and more than one distinct assembly. The local 51 GiB free space cannot hold all shortlist FASTAs uncompressed at once. A full-tree run is therefore a batch/HPC job with sequential cleanup or larger storage, not a local one-shot experiment. Absolute wall time is not defensible without a pilot on the chosen host. A useful extrapolation is `T_full ~= T_pilot * 194.53 Gbp / pilot_Gbp`, adjusted for concurrency. The Takifugu rubripes shortlist is a practical pilot: three accessions but only two distinct assembly keys, totaling 0.776 Gbp. The GCA/GCF Reference mirror may share one result only after FASTA checksum identity is confirmed; EvoSubster still needs exact-accession TSV rows.

## changed_files

- `.conductor/handoffs/T002-handoff.md`

## machine_checks

- PASS: `git status --short --branch` confirmed branch `feature/stage0-reference-quality-ranking`; only pre-existing untracked conductor/user artifacts are visible.
- PASS: executable/version probes confirmed Python 3.9.6, Homebrew 6.0.12, and Java 24.0.2.
- PASS: `command -v` probes confirmed BUSCO and all audited runtime/dependency managers are absent.
- PASS: official dataset metadata probe confirmed `actinopterygii_odb12.2` date, MD5, and archive content length.
- PASS: Stage 0 schema/test and both reported shortlist artifact paths exist and were read successfully.

claim_boundary_ok: true

next_action: Provision the official container on a Linux/Docker host, run the two-distinct-assembly Takifugu pilot first, record timing and checksums, then let the conductor decide whether resource evidence justifies expanding to all 92 potentially overturnable species.

delegation_effect: difficulty=medium; agents=1; useful=identified the missing runtime, exact current lineage artifact, feasible pilot, and full-tree storage boundary; misses=no empirical BUSCO timing because installation/execution was prohibited; next_mix=one execution worker plus one independent TSV/result verifier.
