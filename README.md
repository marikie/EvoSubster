# EvoSubster: Let's investigate evolutionary substitution trends across diverse taxonomic groups ─=≡Σ((( つ•̀ω•́)つ !

## Introduction

We aim to observe single-base and double-base substitution trends across diverse organisms, taking into account the influence of neighboring bases. 

Our pipeline’s input should be NCBI genome accession IDs of three closely related species of your choice (we
suggest >80% identity of orthologous DNA): species A as an outgroup, species B, and species C. It downloads
the corresponding genomic FASTA files and, when available, gene annotations.

Pairwise alignments will be performed between _Species A_ and _Species B_, and between _Species A_ and _Species C_. These two sets of alignments will then be merged into a multiple sequence alignment. Our pipeline infers substitutions in _Species B_ and _Species C_ based on the principle of parsimony, and outputs visualizations.

## Prerequisites

Install the following command-line tools before running any scripts:
- #### NCBI Datasets command-line tools ([https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/))

- #### LAST ([https://gitlab.com/mcfrith/last](https://gitlab.com/mcfrith/last))
- #### yq ([https://github.com/mikefarah/yq](https://github.com/mikefarah/yq))

- #### jq ([https://jqlang.org](https://jqlang.org))

- #### R (≥4.0) with the following libraries:
    * stringr
    * RColorBrewer
    * showtext
    * jsonlite
    * curl
    * dplyr
    * ggplot2
    * rlang
    * sysfonts

- #### Python3 (3.8 or later) with standard library modules

## Repository Setup

```bash
git clone https://github.com/marikie/EvoSubster.git
```

The alignment entry points live under `EvoSubster/last/`, analytics helpers under `EvoSubster/analysis/`, and plotting code under `EvoSubster/analysis/R/`.

## Configuration
### `dwl_config.yaml`
Update `EvoSubster/last/dwl_config.yaml` to point at the directory that should store downloaded genomes:

```yaml
paths:
  # Change the path to the directory to store downloaded genomes
  base_genomes: "/absolute/path/to/genomes"
```
### `sbst_config.yaml`
Adjust `EvoSubster/last/sbst_config.yaml` so results flow to your writable location (This can be overridden with the `--out-dir` flag, later. See below for more details.):

```yaml
paths:
  # Change the path to the directory to store outputs
  out_dir: "/absolute/path/to/outputs"
```

## Running the Pipeline
All commands below must be executed from `EvoSubster/last/`.

### To start from downloading genomes:

Run the following script:

```bash
./trisbst_3spc_fromDwl.sh [DATE] [ORG1_ACCESSION] [ORG2_ACCESSION] [ORG3_ACCESSION] [--out-dir /absolute/path/to/outputs]
```

- `DATE` is an arbitrary run label (e.g., `20250131`).
- The `ORG1` should be the outgroup among the three genomes.  
- `ORG1_ACCESSION`, `ORG2_ACCESSION`, and `ORG3_ACCESSION` are the NCBI accession IDs of the three genomes. (e.g. GCA_023078555.1, GCA_900538255.1, GCA_024500015.1)
- `--out-dir` is optional; when omitted the script writes beneath the `paths.out_dir` configured in `sbst_config.yaml`.

The download wrapper ensures:

- Each accession is fetched via `datasets download genome accession` with the includes defined in `dwl_config.yaml`.
- Archives are unpacked, non-essential folders removed, and FASTA/GFF assets moved to the genome directory.
- `genomic.gff` for the outgroup is auto-detected; if missing, downstream scripts receive `NO_GFF_FILE`.
- FASTA paths and (optional) GFF are passed to `trisbst_3spc.sh`.


### If the genomes are already downloaded:

You can also run the pipeline directly from the downloaded genomes by running the following script:

```bash
./trisbst_3spc.sh [DATE] [ORG1_FASTA] [ORG2_FASTA] [ORG3_FASTA] [ORG1_GFF|NO_GFF_FILE] [--out-dir /absolute/path/to/outputs]
```

- `DATE` is an arbitrary run label (e.g., `20250131`).
- `ORG1_FASTA`, `ORG2_FASTA`, and `ORG3_FASTA` are the paths to the reference FASTA files of the three genomes.
- `ORG1_GFF` is the path to the reference GFF file of the outgroup. If the GFF file is unavailable, enter "NO_GFF_FILE".
- `--out-dir` is optional; when omitted the script writes beneath the `paths.out_dir` configured in `sbst_config.yaml`.

## Outputs

Results are written beneath `paths.out_dir` (or the `--out-dir` override) following this structure:

```
<out_dir>/<ORG1short>_<ORG2short>_<ORG3short>/<DATE>/
```

Key files include:

- `*.train`: you can check the rough substitution percent identity between two species calculated by `LAST` (see [LAST-TRAIN](https://gitlab.com/mcfrith/last/-/blob/main/doc/last-train.rst?ref_type=heads) for more details)
- `*_gcContent*.out`: GC content of the whole genome sequence for the ingroup species
- `*_sbstRatio*.out`: the percentage of single-base substitutions without considering neighboring bases for the ingroup species (see `EvoSubster/analysis/subRatio.py` for more details)
- `*_maflinked.tsv`: single-base substitution counts with maf-link filtering, which removes
alignments between non-homologous insertions of homologous transposons (see [MAF-LINKED](https://gitlab.com/mcfrith/last/-/blob/main/doc/maf-linked.rst?ref_type=heads) of `LAST` for more details)
- `*_maflinked_dinuc.tsv`: dinucleotide substitution counts
- `*_maflinked_sbst.pdf`: bar plot of the number of single-base substitutions
- `*_maflinked_ori.pdf`: bar plot of the number of original trinucleotides
- `*_maflinked_norm.pdf`: bar plot of single-base substitutions normalized by the number of original trinucleotides
- `*_maflinked_logRatio.pdf`: the log₂ of the single-base substitution rate relative to
the overall mean substitution rate across all substitution types (y-axis: log₂[(# of substitutions / # of original trinucleotides) ÷ mean(# of
substitutions / # of original trinucleotides) across all substitution types])
- `*_maflinked_dinuc*.pdf`: bar plot of the number of dinucleotide substitutions normalized by the number of original dinucleotides
- `*_ncds*`: files of non-coding regions analysis

Re-running the pipeline preserves existing outputs and skips recomputation where possible.