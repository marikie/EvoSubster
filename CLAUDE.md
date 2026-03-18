# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EvoSubster is a computational pipeline analyzing evolutionary substitution trends across three closely related organisms (one outgroup, two ingroups). It uses LAST for pairwise alignments, merges them into three-way MAF files, and applies Python/R utilities to examine single-base and dinucleotide substitution patterns.

The core biological model: given three genomes where org1 is the outgroup, substitutions are attributed by parsimony — if the outgroup and one ingroup agree, the other ingroup is assumed to have mutated. Both ingroups produce separate TSV outputs. Edge bases of each trinucleotide/dinucleotide window must be identical across all three genomes, and reverse-complement folding is applied so only C and T midpoints are tracked.

## Build and Run Commands

### Full Pipeline (with genome download)
```bash
./scripts/last/trisbst_3spc_fromDwl.sh <DATE> <ACC1> <ACC2> <ACC3> [--out-dir /path/to/output] [--thread N] [--idt-only]
```
- ACC1 = outgroup accession, ACC2/ACC3 = ingroup accessions

### Pipeline from Existing FASTA
```bash
./scripts/last/trisbst_3spc.sh <DATE> <ORG1.fa> <ORG2.fa> <ORG3.fa> <ORG1.gff|NO_GFF_FILE> [--out-dir /path/to/output] [--thread N] [--idt-only]
```
- `--thread N`: Set number of parallel threads (default: 8)
- `--idt-only`: Run only the percent identity check step

### Bulk Runs
Taxon-specific batch scripts live in `scripts/last/bulkRun_*.sh` (e.g., `bulkRun_fungi.sh`, `bulkRun_actinopteri.sh`). Each invokes `trisbst_3spc_fromDwl.sh` for multiple species triplets.

### Regenerate TSV/Plots Only
```bash
./scripts/last/generate_tsv_files.sh <joined.maf> <out_dir>
./scripts/last/generate_graphs.sh <tsv_dir>
```

### Reporting
```bash
python scripts/last/collect_run_summary.py <data_root> [--output summary.json] [--idt-threshold 80]
bash scripts/last/render_report.sh -j summary.json [-o output.docx] [-f word_document]
```

## Testing

Tests must be run from `scripts/analysis/` (fixtures use relative paths to `./test/`).

```bash
cd scripts/analysis

# Run active test modules
python -m unittest test_subRatio test_trisbst_2TSVs test_trisbst_2TSVs_errprb

# Run a single test class
python -m unittest test_trisbst_2TSVs.TestTriUvMuts2TSVs

# Run a single test method
python -m unittest test_subRatio.TestSubRatio.test_count
```

Test fixtures (`.maf` input and `.tsv` expected output) are in `scripts/analysis/test/`. Tests compare generated output against expected TSVs using `diff`. Older test modules (`test_triSbstTSV`, `test_split2twoFiles`) are in `archive/`.

## Architecture

### Core Data Flow

MAF file → `Util.getJoinedAlignmentObj()` → `JoinedAlignment` objects (with `gSeq1`/`gSeq2`/`gSeq3`) → sliding-window counting → TSV output → R visualization

The `Alignment` class handles pairwise MAF entries (used by alignment utilities); `JoinedAlignment` handles 3-way joined MAF entries (used by all substitution-counting scripts). Both are constructed via `fromMAFEntry()` class methods.

### Key Scripts

**Pipeline orchestration (Bash)** — `scripts/last/`:
- `trisbst_3spc.sh` / `trisbst_3spc_fromDwl.sh` — main entry points
- `last_train.sh` — LAST model training and percent identity check
- `one2one.sh` — pairwise alignment: lastdb → last-train → lastal → last-split → maf-linked
- `mafjoin.sh` — sort and join two pairwise MAFs into a three-way MAF
- `generate_tsv_files.sh` / `generate_graphs.sh` — TSV and plot generation
- `collect_run_summary.py` — gather run metadata into JSON for reporting
- `render_report.sh` / `run_report.sh` — R Markdown report rendering

**Substitution analysis (Python)** — `scripts/analysis/`:
- `trisbst_2TSVs.py` — trinucleotide context substitution counting (main, by M.C. Frith)
- `trisbst_2TSVs_errprb.py` — variant that filters by MAF `p`-line error probability
- `disbst_2TSVs.py` — dinucleotide context substitution counting (by M.C. Frith)
- `subRatio.py` — simple substitution percentage without context
- `isParsimonious.py` — computes parsimony/non-parsimony ratio from percent identities
- `Alignment.py` — `Alignment` (pairwise) and `JoinedAlignment` (3-way) MAF data structures
- `Util.py` — MAF parsing, coordinate conversion, overlap detection, reverse-complement helpers

**Visualization (R)** — `scripts/analysis/R/`:
- `sbmut.R` — substitution bar charts
- `logRatioPlot.R` — log-ratio enrichment plots
- `dinucleotide-plot.R` — dinucleotide substitution plots
- `pca.R` — PCA analysis

### Pipeline Flow
1. GC content calculation
2. LAST training (percent identity between all three pairs)
3. Pairwise one-to-one alignments (org1→org2, org1→org3)
4. `maf-linked` filtering (removes isolated alignments, e.g. non-homologous transposon insertions)
5. Three-way MAF joining via `maf-join`
6. Optional CDS removal via `maf-cut-cds-uglier.py` (produces `_ncds` non-coding variants)
7. Python TSV generation (trinuc, dinuc, with/without maf-linked, with/without ncds)
8. R visualization (bar charts, log-ratio plots, dinucleotide plots)

## Configuration

- `scripts/last/dwl_config.yaml` — genome download paths, default output directory (`paths.out_dir`)
- `scripts/last/sbst_config.yaml` — output filename patterns, messages, settings
- Environment overrides: `YQ_BINARY`, `LAST_DIR_OVERRIDE`, `THREAD_NUM_OVERRIDE`, `THREAD_NUM`

## Dependencies

- **LAST** (lastal, lastdb, last-train, last-split, maf-linked, maf-sort, maf-join)
- **NCBI Datasets CLI** + jq
- **yq** (Mike Farah v4+)
- **Python 3.8+** (stdlib only)
- **R 4.0+** with: stringr, RColorBrewer, showtext, jsonlite, curl, dplyr, ggplot2, rlang, sysfonts

## Coding Conventions

- Python: 4-space indent, snake_case, CPython 3.8+, stdlib only (no third-party deps)
- Shell: `#!/bin/bash`, quote variables, `$(...)` over backticks, idempotent (skip if output exists)
- R: Parameters grouped at top, document external deps inline
- Commits: Imperative mood, optionally scoped (e.g., `last: adjust TSV generation`)

## Git

- Repository root is this directory (`scripts/`). Always run git via `git -C /big/mrk/proj/sbst/scripts`
- **Branching**: Work on feature branches, merge to main when complete
- **Commit granularity**: One commit per feature/change
- **Commit messages**: English, imperative mood, optionally scoped (e.g., `last: add thread option`)
- **Push**: Push directly to GitHub (`marikie/EvoSubster`)
- **`.gitignore`**: Large intermediate files (`*.maf`, `*.train`, etc.) are not tracked. When encountering a new type of large intermediate file, ask the user before adding it to `.gitignore`
- Do not modify files in `archive/` (kept for reference only)
