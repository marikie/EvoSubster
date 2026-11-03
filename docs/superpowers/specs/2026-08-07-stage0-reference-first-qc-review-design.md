# Stage 0 Reference-First QC Review Design

Date: 2026-08-07
Status: Approved in conversation; awaiting review of this written specification

## Purpose

Change eukaryotic Stage 0 from automatic BUSCO-based replacement to a reference-first, audit-first policy. EvoSubster should select an NCBI Reference by default, preferring an annotated Reference when available, retain external BUSCO and Merqury evidence for review, and only replace the baseline assembly when the user explicitly enables a strict QC override.

This policy does not claim that an NCBI Reference is the assembly with the highest base accuracy. It provides a reproducible default suited to coding-sequence analysis and separates assembly selection from quality warnings.

## Scope

The change applies to the eukaryotic external-BUSCO path in Stage 0. Existing prokaryotic ANI, CheckM, and Merqury ranking remains unchanged because the requested problem concerns eukaryotic BUSCO interpretation.

The implementation will update:

- `src/select/assembly_selection.R`
- `src/select/fetch_assembly_metadata.py`
- `src/select/trio_selection.R`
- `src/sbst_fromDwl.sh`
- Stage 0 tests and user documentation

Raw BUSCO runs, Merqury execution, and downstream assembly-sensitivity experiments are outside this change.

## Considered Approaches

### 1. Audit-first default with explicit strict override — selected

Select the highest-priority Reference independently of external BUSCO. Attach QC evidence, identify review cases, and keep the baseline unless the user supplies `--allow-qc-override`. This preserves reproducibility and still supports a deliberate, evidence-bounded replacement.

### 2. Remove BUSCO ranking entirely

This is simpler and prevents accidental replacement, but discards the useful Plectropomus-type case where a distinct assembly strongly improves conserved single-copy gene representation.

### 3. Keep automatic ranking and add thresholds

This retains the current workflow but makes a universal threshold carry too much biological meaning. BUSCO duplication can be technical or biological, and a Complete-only threshold cannot distinguish the Larimichthys failure mode.

## Candidate Eligibility and Baseline Selection

Existing hard exclusions remain: non-current, atypical, ANI Failed, unsupported hybrid/alternate-haplotype/unresolved-diploid, and an optional relative-contig-N50 gate.

For eligible eukaryotic candidates, the baseline order becomes:

1. annotated NCBI Reference;
2. unannotated NCBI Reference;
3. annotated RefSeq assembly;
4. other annotated assembly;
5. existing assembly-level, contiguity, gap, primary/haploid, and accession fallback.

Within one tier, existing metadata quality fields break ties.

The baseline is computed before external QC ranking. Supplying `--assembly-qc` therefore does not change the default selection.

## External QC Schema

The external TSV requires exact versioned GCA/GCF `accession` values and BUSCO
or Merqury evidence on every row. Unmatched rows produce a diagnostic, and an
explicit override fails when the input is absent or empty or when no row
matches the fetched metadata. BUSCO rows use the same genome mode, lineage,
and version within a species. Any row containing a BUSCO-specific value must
set `qc_busco_mode=genome`; Merqury-only rows may leave the mode empty.

Add these optional columns:

- `qc_busco_single`
- `qc_busco_duplicated`
- `qc_busco_internal_stop`

Existing columns remain:

- `qc_busco_complete`
- `qc_busco_fragmented`
- `qc_busco_missing`
- `merqury_qv`
- `merqury_completeness`

BUSCO component fields are percentages in the range 0–100. When Single-copy is absent but Complete and Duplicated are present, Stage 0 derives Single-copy as `Complete - Duplicated`. Stage 0 validates component consistency within normal one-decimal BUSCO rounding tolerance and never compares BUSCO results across species.

Incomplete legacy TSVs are accepted for audit, but cannot trigger an automatic override. The audit explains which components are missing.

## Paired GCA/GCF Assemblies

Add `paired_accession` to fetched metadata when NCBI supplies it. Stage 0 derives an assembly-equivalence key from an accession and its paired accession. A change between paired GCA/GCF accessions is not considered a distinct alternative, does not create an override, and is not reported as a biological winner change.

When paired metadata is unavailable, Stage 0 does not guess equivalence from accession numbers.

## Review Detection

Every candidate receives audit fields:

- `baseline_rank`
- `final_rank`
- `baseline_selected`
- `qc_preferred`
- `review_required`
- `review_reason`
- `override_applied`
- `override_block_reason`
- `assembly_equivalence_key`

Among candidates with complete comparable BUSCO components, `qc_preferred` is determined by Single-copy descending, Duplicated ascending, Missing ascending, Fragmented ascending, followed by the existing assembly metadata order. This comparator proposes a review candidate; it does not change `selected` under the default policy.

A species is marked for review when any of these conditions holds:

- external BUSCO was supplied for the species but the proposed comparison has incomplete BUSCO fields;
- a distinct alternative has stronger comparable single-copy BUSCO evidence;
- the baseline crosses advisory BUSCO warnings: Complete below 90% or Duplicated above 5%;
- an alternative appears better by Complete but is duplication-confounded;
- only one side of a proposed comparison has external QC.

Missing annotation is neither an eligibility exclusion nor a review trigger. The pipeline continues whole-genome analysis, and GFF-dependent non-coding analysis occurs only when a GFF is available.

The 90% and 5% values are review triggers, not exclusion or high-quality guarantees. Absence of external QC alone does not trigger review. Review reasons are semicolon-separated stable codes, with human-readable details retained in the audit columns. The pipeline continues with the baseline selection while review is pending.

`rank_assembly_candidates()` returns a `review` table in addition to `audit`, `shortlist`, and `best`. Tree mode writes `assembly_review.tsv` containing all candidates for species that require review.

## Explicit Override

Add the CLI flag `--allow-qc-override`, default off. `sbst_fromDwl.sh` forwards it to `trio_selection.R`.

With the flag off, the baseline remains selected regardless of external BUSCO or Merqury values.

With the flag on, a distinct alternative may replace the baseline only when all conditions hold:

1. both candidates have comparable BUSCO genome-mode lineage and version;
2. both have Complete, Single-copy or derivable Single-copy, Duplicated, Fragmented, and Missing values;
3. the alternative has strictly higher Single-copy completeness;
4. Duplicated, Fragmented, and Missing are each no worse than the baseline;
5. the candidates do not share an assembly-equivalence key.

Annotation status is not an override precondition; an unannotated alternative may replace the baseline when all conditions above hold.

Merqury values remain supporting audit evidence and are never treated as comparable unless provided for both candidates. They do not rescue a candidate that fails the strict BUSCO dominance conditions. A BUSCO override does not imply verified base-level QV.

If any condition fails, the baseline remains selected and
`override_block_reason` records stable codes for why the proposed override was
not applied, including higher Duplicated, Fragmented, or Missing components.

## Logging and Documentation

Replace the current message that asks users to provide BUSCO for final ranking. Log the number of species requiring review and point to `assembly_review.tsv`. Help and README text must state that external QC is audit-only by default and that `--allow-qc-override` enables strict, explicit replacement.

## Error Handling

- Reject non-genome BUSCO mode.
- Reject unversioned external assembly accessions and diagnose unmatched rows.
- Reject mixed BUSCO lineage or version within a species comparison.
- Reject impossible numeric BUSCO values, including negative values, Complete below Duplicated, or components whose scale is internally inconsistent.
- Accept absent optional QC columns by filling them with missing values and recording incomplete comparison status.
- Never drop the baseline solely because external QC is missing.

## Test Strategy

Tests will be written before production changes and must demonstrate:

1. external BUSCO does not replace an annotated Reference by default;
2. a stronger distinct candidate creates review output;
3. an unannotated Reference remains ahead of an annotated RefSeq assembly;
4. Larimichthys-like increased duplication blocks explicit override;
5. Takifugu-like paired GCA/GCF accessions do not count as an override;
6. Plectropomus-like strict dominance can override when explicitly enabled;
7. incomplete legacy BUSCO fields never override;
8. a species with no annotated candidate remains eligible and creates no annotation-only review warning;
9. the new CLI flag parses and is forwarded by the shell wrapper;
10. existing Stage 0 gate and trio-selection tests remain green.

## Success Criteria

The implementation is complete when the default result is reference-first and annotation-aware, external BUSCO cannot silently change it, all proposed replacements are auditable, paired accessions are not misreported, explicit override is conservative, documentation matches behavior, and all relevant tests pass.
