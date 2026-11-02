# Frozen regression baseline: 20261031101732

These eight expected TSVs are exact copies from the completed pipeline run
`tmp/codex-checks/pipeline-bf779ef6/results/Oikdio1_Oikalb2_Oikvan3/20261031101732/statistics/`.
The user explicitly approved adopting this run as the new baseline. Expectations
must never be regenerated or adjusted automatically from later validation runs.

The corresponding three FASTAs and outgroup GFF are independent, byte-identical
copies under ignored `test_genomes/regression-20261031/`. They are not symlinks to
temporary run data. `manifest.json` records repository-relative input paths,
accessions, SHA256 hashes, FASTA record order, expected TSV hashes, and the source
run's LAST 1648 settings. Source files and the historical 20260420 baseline remain
unchanged. The new baseline preserves the source run's FASTA order: reordering
records can change ambiguous duplicate mappings, including a 279-base joined
block that explains the differences from the historical counts.

Before full validation, activate the repository environment:

```bash
rtk proxy bash -c 'source .venv/bin/activate && python3 scripts/codex_checks.py --queue-pipeline'
```

The worker verifies every input's bytes and FASTA ID order, and checks all eight
expected TSV hashes against the manifest before launching
`sbst.sh` directly. Missing or changed inputs are failures; there is no download
or repair fallback. Restore an exact copy from the approved source snapshot if
verification fails; never sort, rewrite, or replace inputs with fresh downloads.

Input verification, alignment, and counting regression use these frozen local
files without genome downloads or network access. Full figure generation
currently requires Google Fonts network access because `sbmut.R` and
`logRatioPlot.R` call `font_add_google`. A successful counting comparison does
not establish successful figure generation; inspect plotting logs separately.

Full runs explicitly use outgroup annotation, `NO_GFF_FILE` for both ingroups,
the three versioned accessions, eight threads, and `--lastdb-options '-uRY4'`.
Every run uses fresh result, database, and training-cache directories. The log
records verified inputs, interpreter, tool locations/versions, command, and the
byte-for-byte comparison of all eight TSVs (only date filenames are mapped).
The 600-second limit and failure handling apply as before. A timeout or any TSV
mismatch remains a failure; inspect the retained log and `result.json` rather
than editing this baseline.
