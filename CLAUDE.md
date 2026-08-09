# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EvoSubster (PhD thesis: "Evolutionary Patterns of DNA Base Substitutions: Comparative Insights Across Species") is a pipeline for comparative analysis of context-dependent substitution spectra across diverse eukaryotic lineages. It performs whole-genome pairwise alignment with LAST, merges two pairwise MAFs into a three-way MAF, then counts substitutions in trinucleotide context (single-base) and tetranucleotide context (double-base) using parsimony-based Python scripts.

The core biological model: given three genomes where org1 is the outgroup, substitutions are attributed by parsimony — if the outgroup and one ingroup agree, the other ingroup is assumed to have mutated. Both ingroups produce separate TSV outputs. Edge bases of each trinucleotide/dinucleotide window must be identical across all three genomes, and reverse-complement folding is applied so only C and T midpoints are tracked.

Applied to 60 trios across 7 lineages (fungi, cnidarians, Apicomplexa, Phaeophyceae, Oomycota, Porifera, Arthropoda) and reported in the PhD thesis.

**Parsimony case classification:**

- org1 = outgroup, org2/org3 = ingroups (the pipeline's slot convention). The thesis
  and `src/metrics/isParsimonious.py` use A/B/C letters with the outgroup in different
  positions, so do not transcribe their lettering into new code — follow org1/org2/org3.
- Case 0: all three identical → no substitution, increment count
- Case 1: only org1 (outgroup) differs → ambiguous, excluded
- Case 2: only org2 differs → substitution inferred in org2
- Case 3: only org3 differs → substitution inferred in org3
- Case 4: all three differ → multiple scenarios, excluded

## Automation

### Skills (slash commands)
- `/add-trio <lineage> <out_acc> <in1_acc> <in2_acc>` — download genomes + full pipeline
- `/run-lineage <lineage>` — bulk-run all trios in a lineage
- `/gen-figure <lineage> <triplet>` — regenerate R plots from existing TSVs
- `/check-results <lineage> [<triplet>]` — validate output files and identity thresholds
- `/gen-report <lineage> [--format word_document|pdf_document|html_document]` — render summary report

### Hooks (run automatically after file edits)
- `.sh` files: `bash -n` syntax check
- `src/count/`, `src/statistics/`, `src/metrics/` Python files: unit tests
- `src/` or `config/` core files: full pipeline integration test (background, log: `/tmp/pipeline_test_*.log`)

## Repository Layout (Noble 2009)

```
/big/mrk/proj/sbst/evo-subster/        # Project root = Git root
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
├── config/                            # Configuration
│   ├── dwl_config.yaml               # Genome download paths
│   └── sbst_config.yaml              # Output patterns, settings
├── drivers/                           # Bulk execution scripts
│   └── bulkRun_*.sh                   # Per-lineage batch scripts
├── genomes/                           # Reference genome FASTA/GFF files (download target)
├── results/                           # Analysis results (per-lineage)
│   ├── fungi/, cnidaria/, etc.        # Per-lineage output directories
│   └── summary/                       # Cross-lineage aggregation
├── bench_genomes/                     # Benchmark genome inputs
├── bench_results/                     # Benchmark analysis outputs
├── eg_results/                        # Example results for reference
├── tasks/                             # Task tracking (todo.md, lessons.md)
├── doc/                               # Documentation
├── log/                               # Execution logs
├── test_genomes/                      # Integration test genomes (gitignored)
└── test_results/                      # Integration test outputs (gitignored)
```

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

- Repository root is this directory. Git operations can be run directly here.
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
