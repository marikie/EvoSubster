# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EvoSubster (PhD thesis: "Evolutionary Patterns of DNA Base Substitutions: Comparative Insights Across Species") is a pipeline for comparative analysis of context-dependent substitution spectra across diverse eukaryotic lineages. It performs whole-genome pairwise alignment with LAST, merges two pairwise MAFs into a three-way MAF, then counts substitutions in trinucleotide context (single-base) and tetranucleotide context (double-base) using parsimony-based Python scripts.

The core biological model: given three genomes where org1 is the outgroup, substitutions are attributed by parsimony — if the outgroup and one ingroup agree, the other ingroup is assumed to have mutated. Both ingroups produce separate TSV outputs. Edge bases of each trinucleotide/dinucleotide window must be identical across all three genomes, and reverse-complement folding is applied so only C and T midpoints are tracked.

Applied to 60 trios across 7 lineages (fungi, cnidarians, Apicomplexa, Phaeophyceae, Oomycota, Porifera, Arthropoda) and reported in the PhD thesis.

**Parsimony case classification:**

- Species A = outgroup, B/C = ingroups
- Case 0: all three identical → no substitution, increment count
- Case 1: only A differs → ambiguous, excluded
- Case 2: only B differs → substitution inferred in B
- Case 3: only C differs → substitution inferred in C
- Case 4: all three differ → multiple scenarios, excluded

## Repository Layout (Noble 2009)

```
/big/mrk/proj/sbst/                    # Project root = Git root
├── src/                               # Source code (workflow-stage organized)
│   ├── align/                         # LAST alignment + MAF processing
│   │   ├── one2one.sh                 # Pairwise alignment pipeline
│   │   ├── last_train.sh              # LAST model training + identity check
│   │   ├── mafjoin.sh                 # Sort and join two MAFs into three-way
│   │   └── maf-cut-cds-uglier.py      # CDS region filtering
│   ├── count/                         # Substitution counting (core)
│   │   ├── Alignment.py               # Alignment/JoinedAlignment data structures
│   │   ├── Util.py                    # MAF parsing, coordinate conversion
│   │   ├── single_sbst_2TSVs.py       # Trinucleotide SBS counting
│   │   ├── single_sbst_2TSVs_errprob.py # Error probability filtered variant
│   │   └── disbst_2TSVs.py            # Dinucleotide DBS counting
│   ├── statistics/                    # Statistical tests
│   │   ├── chi2_context.py            # Chi-squared context test
│   │   └── chi2_spectrum_cmp.py       # Chi-squared spectrum comparison
│   ├── metrics/                       # Descriptive metrics
│   │   ├── gc_content.sh              # GC content calculation
│   │   ├── subRatio.py                # Substitution ratio
│   │   └── isParsimonious.py          # Parsimony ratio
│   ├── visualize/                     # R visualization scripts
│   │   ├── sbmut.R                    # Substitution bar charts
│   │   ├── logRatioPlot.R             # Log-ratio enrichment plots
│   │   ├── dinucleotide-plot.R        # Dinucleotide plots
│   │   ├── pca.R                      # PCA analysis
│   │   └── pca_2d.R                   # 2D PCA plots
│   ├── report/                        # Report generation
│   │   ├── collect_run_summary.py     # Gather run metadata into JSON
│   │   ├── render_report.sh           # R Markdown report rendering
│   │   ├── run_report.sh              # Collect + render wrapper
│   │   ├── link_tsvs.py              # TSV linking utility
│   │   └── report_template.Rmd       # Report template
│   ├── sbst.sh                        # Main pipeline driver
│   ├── sbst_fromDwl.sh               # Pipeline with genome download
│   ├── dwl_organism.sh               # Genome download utility
│   ├── generate_tsv_files.sh         # TSV regeneration
│   └── generate_graphs.sh            # Graph regeneration
├── test/                              # Tests
│   ├── test_*.py                      # Test modules
│   └── fixtures/                      # Test fixtures (.maf, .tsv)
├── bin/                               # Compiled binaries
├── config/                            # Configuration
│   ├── dwl_config.yaml               # Genome download paths
│   └── sbst_config.yaml              # Output patterns, settings
├── drivers/                           # Bulk execution scripts
│   └── bulkRun_*.sh                   # Per-lineage batch scripts
├── genomes/                           # Reference genome FASTA/GFF files (download target)
├── results/                           # Analysis results (per-lineage)
│   ├── fungi/, cnidaria/, etc.        # Per-lineage output directories
│   └── summary/                       # Cross-lineage aggregation
├── doc/                               # Documentation
├── log/                               # Execution logs
```

## Build and Run Commands

### Full Pipeline (with genome download)

```bash
./src/sbst_fromDwl.sh <DATE> <ACC1> <ACC2> <ACC3> [--out-dir PATH] [--thread N] [--idt-only]
```

- ACC1 = outgroup accession, ACC2/ACC3 = ingroup accessions

### Pipeline from Existing FASTA

```bash
./src/sbst.sh <DATE> <ORG1.fa> <ORG2.fa> <ORG3.fa> <ORG1.gff|NO_GFF_FILE> [--out-dir PATH] [--thread N] [--idt-only]
```

### Bulk Runs

```bash
bash drivers/bulkRun_fungi.sh <DATE>
```

### Regenerate TSV/Plots Only

```bash
./src/generate_tsv_files.sh <joined.maf> <out_dir>
./src/generate_graphs.sh <tsv_dir>
```

### Reporting

```bash
python src/report/collect_run_summary.py <results_root> [--output summary.json] [--idt-threshold 80]
bash src/report/render_report.sh -j summary.json [-o output.docx] [-f word_document]
```

## Testing

Tests must be run from `test/` (fixtures use relative paths to `./fixtures/`).

```bash
cd test

# Run active test modules
python -m unittest test_subRatio test_single_sbst_2TSVs test_single_sbst_2TSVs_errprob

# Run a single test class
python -m unittest test_single_sbst_2TSVs.TestTriUvMuts2TSVs

# Run a single test method
python -m unittest test_subRatio.TestSubRatio.test_count
```

Test fixtures (`.maf` input and `.tsv` expected output) are in `test/fixtures/`. Tests compare generated output against expected TSVs using `diff`.

### Pipeline Verification After Script Changes

After modifying any script under `src/` or `config/`, run the full pipeline (download + analysis) with this test trio to verify nothing is broken:

```bash
rm -rf ./test_genomes
rm -rf ./test_results
now="$(date '+%Y-%m-%d-%H-%M')"
./src/sbst_fromDwl.sh <DATE> GCA_907165135.1 GCA_004367875.1 GCA_004367855.1 --genome-dir ./test_genomes --out-dir ./test_results
```

- `GCA_907165135.1` = outgroup
- `GCA_004367875.1`, `GCA_004367855.1` = ingroups

If the pipeline produces any errors, fix them and re-run until it completes without errors.

## Architecture

### Core Data Flow

MAF file → `Util.getJoinedAlignmentObj()` → `JoinedAlignment` objects (with `gSeq1`/`gSeq2`/`gSeq3`) → sliding-window counting → TSV output → R visualization

### Pipeline Flow

1. GC content calculation (`src/metrics/gc_content.sh`)
2. LAST training — percent identity between all three pairs (`src/align/last_train.sh`)
3. Pairwise one-to-one alignments (`src/align/one2one.sh`)
4. `maf-linked` filtering (removes isolated alignments)
5. Three-way MAF joining via `maf-join` (`src/align/mafjoin.sh`)
6. Optional CDS removal (`src/align/maf-cut-cds-uglier.py`)
7. Python TSV generation (`src/count/`)
8. R visualization (`src/visualize/`)

## Configuration

- `config/dwl_config.yaml` — genome download paths, default output directory (`paths.out_dir`)
- `config/sbst_config.yaml` — output filename patterns, messages, settings
- Environment overrides: `YQ_BINARY`, `THREAD_NUM_OVERRIDE`, `THREAD_NUM`

## Dependencies

- **LAST** (lastal, lastdb, last-train, last-split, maf-linked, maf-sort, maf-join)
- **curl** (for NCBI Datasets API access)
- **jq** (for JSON parsing)
- **yq** (Mike Farah v4+)
- **Python 3.8+** (stdlib only)
- **R 4.0+** with: stringr, RColorBrewer, showtext, jsonlite, curl, dplyr, ggplot2, rlang, sysfonts

## Coding Conventions

- Python: 4-space indent, snake_case, CPython 3.8+, stdlib only (no third-party deps)
- Shell: `#!/bin/bash`, quote variables, `$(...)` over backticks, idempotent (skip if output exists)
- R: Parameters grouped at top, document external deps inline
- Commits: Imperative mood, optionally scoped (e.g., `count: fix trinuc edge case`)

## Git

- Repository root is this directory (project root). Git operations can be run directly here.
- **Branching**: Work on feature branches, merge to main when complete
- **Commit granularity**: One commit per feature/change
- **Commit messages**: English, imperative mood, optionally scoped
- **Push**: Push directly to GitHub (`marikie/EvoSubster`)
- **`.gitignore`**: Large intermediate files (`*.maf`, `*.train`, etc.) are not tracked

## Data File Conventions

- `results/<lineage>/<triplet>/<DATE>/` — analysis output per trio run
- Species are abbreviated (e.g., `oikDio` = _Oikopleura dioica_)
- `.maf` files: pairwise (`one2one`, `many2one`) and joined three-way
- `.train` files: LAST training parameters
- `.tsv` files: substitution count tables
