# Preliminary CDS Annotation and Strand Consistency Check

This report records a preliminary, pre-visualization quality-control analysis of
three-way alignments. The historical filename is retained for task compatibility,
but the analysis does **not** establish orthology. It asks whether regions called
coding in the outgroup alignment row are also called coding in both ingroups and
whether the independently annotated CDS strands are consistent with the relative
MAF row orientations. These results are exploratory and are not a manuscript result.

## 1. Scope

A fresh NCBI Datasets and exact-FTP check on 2026-07-15 found three trios for
which the exact assemblies used by the package have CDS GFFs for both ingroups.
The original task note omitted `puffer` because it listed its B assembly as
`GCA_003711565.1`; the package and MAF use annotated RefSeq assembly
`GCF_003711565.1`.

| scope | trios | reason |
|---|---|---|
| Three-way QC | `lenEdo`, `puffer`, `styPis` | Exact-assembly GFF is available for B and C |
| C-only annotation | `arabidopsis`, `neurospora`, `acrDig` | No exact-assembly B GFF; cannot perform the three-way check |
| No three-way annotation | `agaBis`, `acrMil`, `drosophila`, `pleOst`, `armBor` | At least one, usually both, ingroup GFFs are absent |

The complete availability audit is in
`EvoSubster/config/ncbi_annotation_availability.tsv`. Input URLs, byte counts,
SHA-256 values, and annotation sources are in
`EvoSubster/config/cds_orth_strand_inputs.tsv`.

## 2. Method

### 2.1 Coordinates and CDS segments

GFF CDS rows are converted from one-based inclusive coordinates to zero-based
half-open coordinates. Same-strand overlapping or touching CDS records are
collapsed into disjoint genomic segments, preventing duplicate transcript
isoforms from multiplying the same aligned evidence. Regions covered by both
`+` and `-` CDS records are represented with ambiguous strand `.`.

MAF rows are converted to forward-strand genomic coordinates. For the `k`-th
non-gap base, a `+` row uses `start + k`; a `-` row uses
`seqlen - 1 - (start + k)`. The implementation validates declared MAF spans,
row lengths, strand values, and source-coordinate bounds.

### 2.2 Run definition

A run is a maximal consecutive set of alignment columns in one MAF block where
the A row is non-gap and remains in the same disjoint A CDS segment. An A gap,
non-CDS coordinate, or segment-state change ends the run. The analysis streams
one MAF block and one output row at a time.

### 2.3 Classification

The primary analysis requires at least 10 non-gap bases in both B and C. Runs
below that support, A-strand ambiguity, B/C opposite-strand overlap, or a tied
B/C strand vote are `insufficientEvidence`.

For supported runs, B and C CDS coverage is the fraction of non-gap run bases
covered by any CDS in that species. If either fraction is below 0.5, the run is
`nonCodingIncluded`. Otherwise, each species' CDS strand is compared with its
own MAF row strand. Equal match booleans across A/B/C, whether all true or all
false, are `codingStrandConsistent`; mixed booleans are
`strandContradictions`. All-false is expected when an A minus-strand CDS has the
same relative transcription orientation across the aligned genomes.

## 3. Primary Results

| trio | nonCodingIncluded | codingStrandConsistent | strandContradictions | insufficientEvidence | total runs |
|---|---:|---:|---:|---:|---:|
| lenEdo | 10,283 (14.49%) | 57,972 (81.71%) | 98 (0.14%) | 2,592 (3.65%) | 70,945 |
| puffer | 42,368 (19.12%) | 175,159 (79.03%) | 25 (0.01%) | 4,091 (1.85%) | 221,643 |
| styPis | 15,547 (9.91%) | 135,927 (86.67%) | 44 (0.03%) | 5,314 (3.39%) | 156,832 |

Across the three trios, 449,420 A-CDS runs were emitted. The dominant outcome
was coding/strand consistency. The non-coding-inclusion fraction varied more
substantially among trios than the strand-contradiction fraction. This pattern
does not distinguish biological coding-sequence turnover from annotation or
alignment differences.

## 4. Support-Threshold Sensitivity

| trio | minimum B/C non-gap bases | nonCodingIncluded | codingStrandConsistent | strandContradictions | insufficientEvidence |
|---|---:|---:|---:|---:|---:|
| lenEdo | 1 | 11,852 | 58,921 | 113 | 59 |
| lenEdo | 10 | 10,283 | 57,972 | 98 | 2,592 |
| lenEdo | 30 | 7,914 | 54,129 | 90 | 8,812 |
| puffer | 1 | 44,427 | 177,111 | 27 | 78 |
| puffer | 10 | 42,368 | 175,159 | 25 | 4,091 |
| puffer | 30 | 38,529 | 171,086 | 22 | 12,006 |
| styPis | 1 | 18,294 | 138,415 | 47 | 76 |
| styPis | 10 | 15,547 | 135,927 | 44 | 5,314 |
| styPis | 30 | 11,888 | 131,329 | 36 | 13,579 |

The exact values are in `data/results/tables/cds_orth_strand_sensitivity.csv`.
The primary threshold of 10 excludes short alignments without making the
headline counts depend on a single-base majority vote.

## 5. Concrete Examples

### 5.1 `codingStrandConsistent`

| trio | block | run | A segment/strand | B covered/non-gap, strand | C covered/non-gap, strand |
|---|---:|---:|---|---|---|
| lenEdo | 1 | 0 | `NW_025764307.1:4854-7031 (-)` | 48/48 (`+`; MAF `-`) | 48/48 (`+`; MAF `-`) |
| puffer | 39 | 1 | `NC_018890.1:230073-230264 (+)` | 191/191 (`-`; MAF `-`) | 191/191 (`+`; MAF `+`) |
| styPis | 6 | 0 | `NW_019217784.1:7236-7359 (+)` | 123/123 (`-`; MAF `-`) | 123/123 (`-`; MAF `-`) |

### 5.2 `nonCodingIncluded`

| trio | block | run | A segment/strand | B covered/non-gap | C covered/non-gap |
|---|---:|---:|---|---:|---:|
| lenEdo | 3 | 2 | `NW_025764307.1:7923-8119 (-)` | 43/196 | 196/196 |
| puffer | 37 | 0 | `NC_018890.1:218682-218842 (+)` | 160/160 | 0/160 |
| styPis | 0 | 0 | `NC_011162.1:2007-2718 (+)` | 0/711 | 711/711 |

### 5.3 `strandContradictions`

| trio | block | run | A strand/MAF | B strand/MAF | C strand/MAF |
|---|---:|---:|---|---|---|
| lenEdo | 79 | 0 | `+ / +` | `+ / -` | `- / +` |
| puffer | 655 | 41 | `- / +` | `+ / -` | `- / -` |
| styPis | 8,516 | 0 | `+ / +` | `- / +` | `+ / -` |

### 5.4 `insufficientEvidence`

| trio | block | run | B covered/non-gap | C covered/non-gap | reason |
|---|---:|---:|---:|---:|---|
| lenEdo | 3 | 1 | 19/19 | 0/7 | B or C has fewer than 10 non-gap bases |
| puffer | 41 | 3 | 0/3 | 0/3 | B or C has fewer than 10 non-gap bases |
| styPis | 3 | 1 | 0/3 | 3/3 | B or C has fewer than 10 non-gap bases |

## 6. Verification and Performance

- All four class counts sum to `total_runs` for every trio and support threshold.
- Each raw CSV row count and independently recounted class distribution match
  its summary row. No MAF blocks had fewer than three rows.
- Nine rows were independently spot-checked without importing the analysis
  module: one consistent, one non-coding, and one strand-contradiction row per
  trio. Direct MAF coordinate conversion and GFF CDS/strand recounting matched
  every reported value and classification.
- `lenEdo` completed in 57.5 seconds with about 100 MB peak RSS.
- `puffer` completed in 642.6 seconds with about 272 MB peak RSS.
- The final combined run completed in 812.1 seconds with about 270 MB peak RSS.

## 7. Interpretation Limits

`codingStrandConsistent` means that independent CDS calls and relative strand
orientations agree at the aligned locus. It does not establish one-to-one
orthology, conserved synteny, absence of paralogy, or correctness of any
individual annotation. `nonCodingIncluded` combines possible biological
coding-sequence differences with annotation incompleteness and alignment
projection errors. These categories should guide inspection before
visualization, not be promoted directly as biological conclusions.

## 8. File Map

- `data/alignments/{lenEdo,puffer,styPis}_ingroup{B,C}_CDS.gff.gz`: six exact-assembly independent ingroup annotations.
- `data/results/tables/cds_orth_strand_raw_<trio>.csv`: run-level evidence for each trio.
- `data/results/tables/cds_orth_strand_summary.csv`: primary threshold-10 counts.
- `data/results/tables/cds_orth_strand_sensitivity.csv`: threshold 1/10/30 counts.
- `EvoSubster/src/analysis/cds_orth_strand_check.py`: streaming classifier and CLI.
- `EvoSubster/src/analysis/cds_intervals.py`: disjoint CDS strand-state index.
- `EvoSubster/src/align/maf_io.py`: strict MAF parser and forward-coordinate helpers.
- `EvoSubster/config/cds_orth_strand_inputs.tsv`: input provenance and checksums.
- `EvoSubster/config/ncbi_annotation_availability.tsv`: 11-trio scope audit.
- `EvoSubster/test/test_cds_orth_strand_check.py`, `test_cds_intervals.py`, and `test_maf_io.py`: synthetic and unit tests.

## 9. Reproduction

From the `EvoSubster` repository root:

```bash
python3 src/analysis/cds_orth_strand_check.py \
  --data-root "/absolute/path/to/Codon-aware coding spectra package" \
  --trios lenEdo puffer styPis \
  --min-nongap 10 \
  --sensitivity-thresholds 1,10,30
```

The CLI verifies the versioned manifest, SHA-256 values, assembly accession
directives, CDS presence, and MAF/GFF seqid compatibility before classification.
