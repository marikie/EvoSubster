#!/usr/bin/env Rscript

# Behavior tests for the two-phase Stage 0 assembly selection:
#   1. hard quality exclusions and reference-anchored shortlisting;
#   2. prokaryote CheckM/ANI ranking or eukaryote reference-first baseline ranking.
# Runs entirely on synthetic metadata and never contacts NCBI.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))
selection_module <- file.path(root, "src", "select", "assembly_selection.R")
if (!file.exists(selection_module)) {
  cat("FAIL: Stage 0 assembly-selection module does not exist\n")
  quit(status = 1)
}
source(selection_module)

fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) cat("ok:", msg, "\n") else { cat("FAIL:", msg, "\n"); fail <<- fail + 1L }
}

row <- function(accession, species, reference = FALSE, level = "Chromosome",
                status = "current", atypical = FALSE, assembly_type = "haploid",
                atypical_warning = "",
                ani = "", checkm_complete = NA_real_, checkm_contam = NA_real_,
                ncbi_busco_complete = NA_real_, contigs = 10, contig_n50 = 1e6,
                total_length = 1e8, ungapped_length = 0.99e8,
                source_database = if (reference) "SOURCE_DATABASE_REFSEQ" else "SOURCE_DATABASE_GENBANK",
                annotated = reference, paired_accession = "") {
  data.frame(
    accession = accession,
    current_accession = accession,
    organism_name = species,
    organism_tax_id = if (grepl("^Prok", species)) "100" else "200",
    species = species,
    genus = sub(" .*", "", species),
    source_database = source_database,
    refseq_category = if (reference) "reference genome" else "",
    assembly_level = level,
    assembly_type = assembly_type,
    release_date = "2025-01-01",
    sequencing_tech = "PacBio",
    assembly_method = "ExampleAssembler",
    is_atypical = if (atypical) "true" else "false",
    atypical_warnings = if (nzchar(atypical_warning)) atypical_warning else
      if (atypical) "contaminated" else "",
    ani_check_status = ani,
    checkm_completeness = checkm_complete,
    checkm_contamination = checkm_contam,
    busco_lineage = if (is.na(ncbi_busco_complete)) "" else "euk_odb10",
    busco_version = if (is.na(ncbi_busco_complete)) "" else "5.7.1",
    busco_complete = ncbi_busco_complete,
    busco_duplicated = NA_real_,
    busco_fragmented = NA_real_,
    busco_missing = NA_real_,
    total_sequence_length = total_length,
    total_ungapped_length = ungapped_length,
    number_of_contigs = contigs,
    contig_n50 = contig_n50,
    number_of_scaffolds = contigs,
    scaffold_n50 = contig_n50,
    has_annotation = if (annotated) "true" else "false",
    annotation_provider = if (annotated) "NCBI RefSeq" else "",
    annotation_release_date = if (annotated) "2025-02-01" else "",
    paired_accession = paired_accession,
    assembly_status = status,
    stringsAsFactors = FALSE
  )
}

meta <- do.call(rbind, list(
  row("P_ref", "Prok alpha", reference = TRUE, level = "Complete Genome",
      ani = "OK", checkm_complete = 98.0, checkm_contam = 0.5, contigs = 1,
      contig_n50 = 4.8e6, total_length = 4.8e6, ungapped_length = 4.8e6),
  row("P_best", "Prok alpha", ani = "OK", checkm_complete = 99.5,
      checkm_contam = 0.1, contigs = 2, contig_n50 = 3e6,
      total_length = 4.9e6, ungapped_length = 4.9e6),
  row("P_ani_failed", "Prok alpha", ani = "Failed", checkm_complete = 100,
      checkm_contam = 0),
  row("P_atypical", "Prok alpha", atypical = TRUE, ani = "OK",
      checkm_complete = 100, checkm_contam = 0),
  row("P_single_cell", "Prok alpha", atypical_warning = "derived from single cell",
      ani = "OK", checkm_complete = 100, checkm_contam = 0),
  row("P_old", "Prok alpha", status = "suppressed", ani = "OK",
      checkm_complete = 100, checkm_contam = 0),
  row("P2_A_low_qv", "Prok delta", ani = "OK", checkm_complete = 99,
      checkm_contam = 0.1, contigs = 2, contig_n50 = 3e6),
  row("P2_Z_high_qv", "Prok delta", ani = "OK", checkm_complete = 99,
      checkm_contam = 0.1, contigs = 2, contig_n50 = 3e6),
  row("E_ref", "Euk beta", reference = TRUE, ncbi_busco_complete = 0.999,
      contig_n50 = 2e6),
  row("E_qc", "Euk beta", ncbi_busco_complete = 0.90, contig_n50 = 5e6),
  row("E_unresolved", "Euk beta", assembly_type = "unresolved diploid",
      ncbi_busco_complete = 1.0, contig_n50 = 10e6),
  row("E2_ref", "Euk gamma", reference = TRUE, ncbi_busco_complete = 0.90,
      contig_n50 = 1e6),
  row("E2_long", "Euk gamma", ncbi_busco_complete = 0.99, contig_n50 = 8e6),
  row("E3_reference", "Euk delta", reference = TRUE, annotated = FALSE,
      ncbi_busco_complete = 0.80, contig_n50 = 1e6),
  row("GCF_003_refseq", "Euk delta", annotated = TRUE,
      source_database = "SOURCE_DATABASE_GENBANK", ncbi_busco_complete = 0.99,
      contig_n50 = 8e6),
  row("RS_004_refseq", "Euk epsilon", annotated = TRUE,
      source_database = "SOURCE_DATABASE_REFSEQ", ncbi_busco_complete = 0.85,
      contig_n50 = 1e6),
  row("GCA_004_other", "Euk epsilon", annotated = TRUE,
      ncbi_busco_complete = 0.99, contig_n50 = 8e6),
  row("E6_unannotated", "Euk zeta", annotated = FALSE,
      ncbi_busco_complete = 0.95, contig_n50 = 8e6,
      paired_accession = "GCF_006_paired")
))
meta$source_database[meta$accession == "GCA_004_other"] <- NA_character_

external_qc <- data.frame(
  accession = c("E_ref", "E_qc", "P2_A_low_qv", "P2_Z_high_qv"),
  qc_busco_mode = c("genome", "genome", "", ""),
  qc_busco_lineage = c("eukaryota_odb10", "eukaryota_odb10", "", ""),
  qc_busco_version = c("5.8.2", "5.8.2", "", ""),
  qc_busco_complete = c(95, 99, NA, NA),
  qc_busco_fragmented = c(2, 0.5, NA, NA),
  qc_busco_missing = c(3, 0.5, NA, NA),
  merqury_qv = c(45, 50, 30, 60),
  merqury_completeness = c(96, 99, 95, 99),
  stringsAsFactors = FALSE
)

ranked <- rank_assembly_candidates(meta, min_rel_contig_n50 = 0,
                                   shortlist_k = 3, external_qc = external_qc)
best <- ranked$best
pick <- function(sp) best$accession[best$species == sp]

ranked_k1 <- rank_assembly_candidates(meta, min_rel_contig_n50 = 0,
                                      shortlist_k = 1, external_qc = external_qc)
best_k1 <- ranked_k1$best

check(identical(pick("Prok alpha"), "P_best"),
      "prokaryote: higher CheckM completeness and lower contamination can beat reference")
check(identical(best_k1$accession[best_k1$species == "Prok alpha"], "P_ref"),
      "prokaryote: shortlist_k=1 selects the phase-1 reference, not an outside candidate")
check(identical(pick("Prok delta"), "P2_Z_high_qv"),
      "prokaryote: Merqury breaks a tie between otherwise equivalent final candidates")
check(identical(pick("Euk beta"), "E_ref"),
      "eukaryote: external QC does not replace an annotated Reference by default")
check(identical(pick("Euk gamma"), "E2_ref"),
      "eukaryote without external genome QC: reference remains the provisional choice")
check(identical(pick("Euk delta"), "E3_reference"),
      "eukaryote tier: unannotated Reference beats annotated non-reference RefSeq")
check(identical(pick("Euk epsilon"), "RS_004_refseq"),
      "eukaryote tier: annotated RefSeq beats another annotated assembly")

audit <- ranked$audit
reason <- function(acc) audit$exclusion_reason[audit$accession == acc]
check(identical(reason("P_ani_failed"), "ani_failed"),
      "hard filter: ANI Failed candidate is excluded with an audit reason")
check(identical(reason("P_atypical"), "atypical: contaminated"),
      "hard filter: atypical candidate is excluded with its NCBI warning")
check(identical(reason("P_single_cell"), "atypical: derived from single cell"),
      "hard filter: atypical source-material warning is excluded before ranking")
check(identical(reason("P_old"), "not_current"),
      "hard filter: suppressed/replaced candidate is excluded")
check(identical(reason("E_unresolved"), "unsupported_assembly_type: unresolved diploid"),
      "hard filter: unresolved diploid candidate is excluded")

short_prok <- ranked$shortlist[ranked$shortlist$species == "Prok alpha", ]
check("P_ref" %in% short_prok$accession,
      "phase 1: eligible NCBI reference is always retained in the shortlist")
check(nrow(short_prok) <= 3,
      "phase 1: shortlist is capped at the requested top-k per species")

euk_basis <- best$selection_basis[best$species == "Euk beta"]
check(identical(euk_basis, "reference_metadata_baseline"),
      "audit: eukaryote selection records the metadata baseline basis")

required_audit_fields <- c(
  "baseline_selected", "qc_preferred", "review_required", "review_reason",
  "override_applied", "assembly_equivalence_key"
)
check(all(required_audit_fields %in% names(audit)),
      "audit: reference-first selection initializes all review-state fields")
check(isTRUE(audit$is_refseq[audit$accession == "GCF_003_refseq"]) &&
        isTRUE(audit$is_refseq[audit$accession == "RS_004_refseq"]) &&
        identical(audit$is_refseq[audit$accession == "GCA_004_other"], FALSE),
      "audit: RefSeq detection accepts GCF accessions and source metadata without NA values")
check(identical(
  audit$assembly_equivalence_key[audit$accession == "E6_unannotated"],
  "E6_unannotated|GCF_006_paired"
), "audit: paired accessions use a stable lexical equivalence key")
check(identical(audit$selected, audit$baseline_selected),
      "audit: default selection remains identical to the metadata baseline")
check(is.data.frame(ranked$review) &&
        !("Euk zeta" %in% ranked$review$species),
      "review: missing annotation alone does not create a review record")

bad_qc <- external_qc
bad_qc$qc_busco_lineage[2] <- "metazoa_odb10"
err <- try(rank_assembly_candidates(meta, external_qc = bad_qc), silent = TRUE)
check(inherits(err, "try-error") && grepl("same BUSCO lineage", as.character(err)),
      "external QC: mixed BUSCO lineages within one species are rejected")

incomplete_qc <- external_qc
incomplete_qc$qc_busco_version[2] <- ""
err <- try(rank_assembly_candidates(meta, external_qc = incomplete_qc), silent = TRUE)
check(inherits(err, "try-error") && grepl("lineage and version", as.character(err)),
      "external QC: BUSCO genome-mode rows must record lineage and version")

strict_qc_row <- function(accession, complete, single = NA_real_, duplicated = NA_real_,
                          fragmented = NA_real_, missing = NA_real_,
                          internal_stop = NA_real_) {
  data.frame(
    accession = accession,
    qc_busco_mode = "genome",
    qc_busco_lineage = "actinopterygii_odb10",
    qc_busco_version = "5.8.2",
    qc_busco_complete = complete,
    qc_busco_single = single,
    qc_busco_duplicated = duplicated,
    qc_busco_fragmented = fragmented,
    qc_busco_missing = missing,
    qc_busco_internal_stop = internal_stop,
    merqury_qv = NA_real_,
    merqury_completeness = NA_real_,
    stringsAsFactors = FALSE
  )
}

strict_pair <- function(species, baseline_accession, alternative_accession,
                        baseline_paired = "", alternative_paired = "") {
  rbind(
    row(baseline_accession, species, reference = TRUE, annotated = TRUE,
        ncbi_busco_complete = 0.95, paired_accession = baseline_paired),
    row(alternative_accession, species, annotated = FALSE,
        ncbi_busco_complete = 0.96, contig_n50 = 3e6,
        paired_accession = alternative_paired)
  )
}

derived_meta <- strict_pair("Euk derived", "DER_ref", "DER_alt")
derived_qc <- rbind(
  strict_qc_row("DER_ref", 96, duplicated = 2, fragmented = 1, missing = 3),
  strict_qc_row("DER_alt", 97, duplicated = 1, fragmented = 1, missing = 2)
)
derived <- rank_assembly_candidates(derived_meta, external_qc = derived_qc)$audit
check(identical(derived$qc_busco_single, c(94, 96)),
      "external QC: Single-copy is derived as Complete minus Duplicated")

invalid_qc_cases <- list(
  list(
    label = "non-numeric percentage",
    qc = within(
      strict_qc_row("VAL_ref", 95, 93, 2, 2, 3),
      qc_busco_single <- "not-a-number"
    ),
    field = "qc_busco_single"
  ),
  list(
    label = "percentage above 100",
    qc = strict_qc_row("VAL_ref", 101, 99, 2, 0, 0),
    field = "qc_busco_complete"
  ),
  list(
    label = "percentage below 0",
    qc = strict_qc_row("VAL_ref", 99, 100, -1, 0, 1),
    field = "qc_busco_duplicated"
  ),
  list(
    label = "Complete below Duplicated",
    qc = strict_qc_row("VAL_ref", 40, 0, 50, 10, 50),
    field = "qc_busco_complete"
  ),
  list(
    label = "Complete-Single-Duplicated inconsistency",
    qc = strict_qc_row("VAL_ref", 95, 90, 2, 2, 3),
    field = "qc_busco_single"
  ),
  list(
    label = "Complete-Fragmented-Missing inconsistency",
    qc = strict_qc_row("VAL_ref", 95, 93, 2, 1, 1),
    field = "qc_busco_complete"
  ),
  list(
    label = "internal-stop percentage above 100",
    qc = strict_qc_row("VAL_ref", 95, 93, 2, 2, 3, 101),
    field = "qc_busco_internal_stop"
  )
)
validation_meta <- row("VAL_ref", "Euk validation", reference = TRUE,
                       ncbi_busco_complete = 0.95)
for (case in invalid_qc_cases) {
  err <- try(rank_assembly_candidates(validation_meta, external_qc = case$qc),
             silent = TRUE)
  check(inherits(err, "try-error") && grepl("VAL_ref", as.character(err)) &&
          grepl(case$field, as.character(err), fixed = TRUE),
        paste("external QC validation:", case$label,
              "identifies accession and field"))
}

rounding_boundary_qc <- strict_qc_row(
  "VAL_ref", 95, 92.8, 2, 2.2, 3
)
err <- try(rank_assembly_candidates(
  validation_meta, external_qc = rounding_boundary_qc
), silent = TRUE)
check(!inherits(err, "try-error"),
      "external QC validation: an exact 0.2 BUSCO rounding difference is accepted")

legacy_meta <- strict_pair("Euk legacy", "LEG_ref", "LEG_alt")
legacy_qc <- rbind(
  strict_qc_row("LEG_ref", 95, fragmented = 2, missing = 3),
  strict_qc_row("LEG_alt", 98, fragmented = 1, missing = 1)
)
legacy <- rank_assembly_candidates(
  legacy_meta, external_qc = legacy_qc, allow_qc_override = TRUE
)
legacy_reason <- paste(legacy$audit$review_reason, collapse = ";")
check(identical(legacy$best$accession, "LEG_ref") &&
        grepl("incomplete_busco_comparison", legacy_reason, fixed = TRUE),
      "override: incomplete legacy BUSCO rows are accepted for audit but cannot override")

plectropomus_meta <- strict_pair(
  "Plectropomus leopardus", "PLE_ref", "PLE_alt"
)
plectropomus_qc <- rbind(
  strict_qc_row("PLE_ref", 96, 94, 2, 1, 3),
  strict_qc_row("PLE_alt", 98, 97, 1, 0.5, 1.5)
)
plectropomus_default <- rank_assembly_candidates(
  plectropomus_meta, external_qc = plectropomus_qc
)
plectropomus_override <- rank_assembly_candidates(
  plectropomus_meta, external_qc = plectropomus_qc,
  allow_qc_override = TRUE
)
plectropomus_audit <- plectropomus_override$audit
check(identical(plectropomus_default$best$accession, "PLE_ref") &&
        any(plectropomus_default$audit$qc_preferred[
          plectropomus_default$audit$accession == "PLE_alt"
        ]),
      "review: a distinct higher-Single-copy candidate is QC-preferred without changing the default")
check(identical(plectropomus_override$best$accession, "PLE_alt") &&
        isTRUE(plectropomus_audit$override_applied[
          plectropomus_audit$accession == "PLE_alt"
        ]) &&
        identical(plectropomus_override$best$selection_basis,
                  "explicit_busco_override"),
      "override: explicit flag permits strict BUSCO dominance by a distinct candidate")
check(any(grepl(
  "alternative_higher_single_copy",
  plectropomus_default$audit$review_reason,
  fixed = TRUE
)), "review: strict higher-Single-copy evidence records a stable reason code")

larimichthys_meta <- strict_pair(
  "Larimichthys crocea", "LAR_ref", "LAR_alt"
)
larimichthys_qc <- rbind(
  strict_qc_row("LAR_ref", 95, 94, 1, 2, 3),
  strict_qc_row("LAR_alt", 97, 93, 4, 1, 2)
)
larimichthys <- rank_assembly_candidates(
  larimichthys_meta, external_qc = larimichthys_qc,
  allow_qc_override = TRUE
)
check(identical(larimichthys$best$accession, "LAR_ref") && any(grepl(
  "complete_gain_duplication_confounded",
  larimichthys$audit$review_reason,
  fixed = TRUE
)), "override: a Complete gain with worse duplication stays baseline and is reviewable")

takifugu_meta <- strict_pair(
  "Takifugu rubripes", "GCF_000180615.1", "GCA_000180615.2",
  baseline_paired = "GCA_000180615.2",
  alternative_paired = "GCF_000180615.1"
)
takifugu_qc <- rbind(
  strict_qc_row("GCF_000180615.1", 96, 95, 1, 1, 3),
  strict_qc_row("GCA_000180615.2", 98, 97, 1, 0.5, 1.5)
)
takifugu <- rank_assembly_candidates(
  takifugu_meta, external_qc = takifugu_qc, allow_qc_override = TRUE
)
check(identical(takifugu$best$accession, "GCF_000180615.1") &&
        !any(takifugu$audit$review_required) &&
        !any(takifugu$audit$override_applied),
      "paired accessions: equivalent GCF/GCA rows neither trigger review nor override")

one_sided_meta <- strict_pair("Euk one-sided", "ONE_ref", "ONE_alt")
one_sided_qc <- strict_qc_row("ONE_alt", 98, 97, 1, 0.5, 1.5)
one_sided <- rank_assembly_candidates(
  one_sided_meta, external_qc = one_sided_qc, allow_qc_override = TRUE
)
check(identical(one_sided$best$accession, "ONE_ref") && any(grepl(
  "one_sided_external_qc", one_sided$audit$review_reason, fixed = TRUE
)), "review: one-sided external QC is recorded and cannot override")

warning_meta <- rbind(
  strict_pair("Euk low-complete", "LOW_ref", "LOW_alt"),
  strict_pair("Euk high-duplicate", "DUP_ref", "DUP_alt")
)
warning_qc <- do.call(rbind, list(
  strict_qc_row("LOW_ref", 89, 87, 2, 3, 8),
  strict_qc_row("LOW_alt", 88, 86, 2, 3, 9),
  strict_qc_row("DUP_ref", 96, 90, 6, 1, 3),
  strict_qc_row("DUP_alt", 95, 90, 5, 1, 4)
))
warning_audit <- rank_assembly_candidates(
  warning_meta, external_qc = warning_qc
)$audit
check(any(grepl("baseline_busco_complete_below_90",
                warning_audit$review_reason[warning_audit$species == "Euk low-complete"],
                fixed = TRUE)),
      "review: baseline Complete below 90 records the advisory warning")
check(any(grepl("baseline_busco_duplicated_above_5",
                warning_audit$review_reason[warning_audit$species == "Euk high-duplicate"],
                fixed = TRUE)),
      "review: baseline Duplicated above 5 records the advisory warning")

single_warning_meta <- row(
  "SOLO_ref", "Euk solo warning", reference = TRUE,
  ncbi_busco_complete = 0.89
)
single_warning_qc <- strict_qc_row("SOLO_ref", 89, 87, 2, 3, 8)
single_warning_audit <- rank_assembly_candidates(
  single_warning_meta, external_qc = single_warning_qc
)$audit
check(grepl(
  "baseline_busco_complete_below_90",
  single_warning_audit$review_reason,
  fixed = TRUE
), "review: baseline BUSCO warnings do not require a distinct alternative")

qc_order_meta <- do.call(rbind, list(
  strict_pair("Euk QC D", "QD_ref", "QD_alt"),
  strict_pair("Euk QC M", "QM_ref", "QM_alt"),
  strict_pair("Euk QC F", "QF_ref", "QF_alt"),
  strict_pair("Euk QC metadata", "QX_ref", "QX_alt")
))
qc_order_qc <- do.call(rbind, list(
  strict_qc_row("QD_ref", 95, 90, 5, 2, 3),
  strict_qc_row("QD_alt", 94, 90, 4, 3, 3),
  strict_qc_row("QM_ref", 94, 90, 4, 2, 4),
  strict_qc_row("QM_alt", 94, 90, 4, 3, 3),
  strict_qc_row("QF_ref", 94, 90, 4, 2, 4),
  strict_qc_row("QF_alt", 94, 90, 4, 1.9, 4),
  strict_qc_row("QX_ref", 95, 92, 3, 2, 3),
  strict_qc_row("QX_alt", 95, 92, 3, 2, 3)
))
qc_order_audit <- rank_assembly_candidates(
  qc_order_meta, shortlist_k = Inf, external_qc = qc_order_qc
)$audit
qc_pick <- function(species) {
  qc_order_audit$accession[
    qc_order_audit$species == species & qc_order_audit$qc_preferred
  ]
}
check(identical(qc_pick("Euk QC D"), "QD_alt"),
      "QC preference: lower Duplicated breaks a Single-copy tie")
check(identical(qc_pick("Euk QC M"), "QM_alt"),
      "QC preference: lower Missing follows Duplicated")
check(identical(qc_pick("Euk QC F"), "QF_alt"),
      "QC preference: lower Fragmented follows Missing")
check(identical(qc_pick("Euk QC metadata"), "QX_ref"),
      "QC preference: existing metadata order is the final tie-breaker")

strict_blocker_meta <- rbind(
  strict_pair("Euk fragmented blocker", "FB_ref", "FB_alt"),
  strict_pair("Euk missing blocker", "MB_ref", "MB_alt")
)
strict_blocker_qc <- do.call(rbind, list(
  strict_qc_row("FB_ref", 95, 94, 1, 2, 3),
  strict_qc_row("FB_alt", 96, 95, 1, 3, 1),
  strict_qc_row("MB_ref", 95, 94, 1, 2, 3),
  strict_qc_row("MB_alt", 96, 95, 1, 0, 4)
))
strict_blockers <- rank_assembly_candidates(
  strict_blocker_meta, external_qc = strict_blocker_qc,
  allow_qc_override = TRUE
)$best
check(identical(
  strict_blockers$accession[
    strict_blockers$species == "Euk fragmented blocker"
  ],
  "FB_ref"
), "override: worse Fragmented blocks a higher-Single-copy alternative")
check(identical(
  strict_blockers$accession[
    strict_blockers$species == "Euk missing blocker"
  ],
  "MB_ref"
), "override: worse Missing blocks a higher-Single-copy alternative")

if (fail > 0L) { cat("\n", fail, " test(s) FAILED\n", sep = ""); quit(status = 1) }
cat("\nAll Stage 0 quality-selection tests passed.\n")
