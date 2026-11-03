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
  row("GCA_900000101.1", "Prok delta", ani = "OK", checkm_complete = 99,
      checkm_contam = 0.1, contigs = 2, contig_n50 = 3e6),
  row("GCA_900000102.1", "Prok delta", ani = "OK", checkm_complete = 99,
      checkm_contam = 0.1, contigs = 2, contig_n50 = 3e6),
  row("GCF_900000201.1", "Euk beta", reference = TRUE, ncbi_busco_complete = 0.999,
      contig_n50 = 2e6),
  row("GCA_900000202.1", "Euk beta", ncbi_busco_complete = 0.90, contig_n50 = 5e6),
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
  accession = c(
    "GCF_900000201.1", "GCA_900000202.1",
    "GCA_900000101.1", "GCA_900000102.1"
  ),
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
check(identical(pick("Prok delta"), "GCA_900000102.1"),
      "prokaryote: Merqury breaks a tie between otherwise equivalent final candidates")
check(identical(pick("Euk beta"), "GCF_900000201.1"),
      "eukaryote: external QC does not replace an annotated Reference by default")
check(identical(pick("Euk gamma"), "E2_ref"),
      "eukaryote without external genome QC: reference remains the metadata baseline")
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
  "baseline_rank", "baseline_selected", "qc_preferred", "review_required",
  "review_reason", "override_applied", "override_block_reason",
  "assembly_equivalence_key"
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

profile_meta <- rbind(
  row(
    "GCA_910000001.1", "Euk profile invariant", annotated = FALSE,
    level = "Complete Genome", contig_n50 = 8e6
  ),
  row(
    "GCF_910000002.1", "Euk profile invariant", annotated = TRUE,
    source_database = "SOURCE_DATABASE_REFSEQ", level = "Contig",
    contig_n50 = 1e5
  )
)
profile_qc <- rbind(
  strict_qc_row("GCA_910000001.1", 95, 94, 1, 2, 3),
  strict_qc_row("GCF_910000002.1", 96, 95, 1, 1, 3)
)
profile_without_qc <- rank_assembly_candidates(
  profile_meta, shortlist_k = Inf
)$audit
profile_with_qc <- rank_assembly_candidates(
  profile_meta, shortlist_k = Inf, external_qc = profile_qc
)$audit
selected_accessions <- function(audit_table, field) {
  sort(audit_table$accession[audit_table[[field]]])
}
check(
  identical(
    selected_accessions(profile_without_qc, "baseline_selected"),
    selected_accessions(profile_with_qc, "baseline_selected")
  ) &&
    identical(
      selected_accessions(profile_without_qc, "selected"),
      selected_accessions(profile_with_qc, "selected")
    ) &&
    identical(
      selected_accessions(profile_without_qc, "selected"),
      "GCF_910000002.1"
    ),
  "policy invariant: external QC cannot change a non-prokaryote baseline or default selection"
)

missing_mode_qc <- strict_qc_row(
  "GCF_910000002.1", 96, 95, 1, 1, 3
)
missing_mode_qc$qc_busco_mode <- ""
err <- try(
  rank_assembly_candidates(profile_meta, external_qc = missing_mode_qc),
  silent = TRUE
)
check(
  inherits(err, "try-error") &&
    grepl("qc_busco_mode=genome", as.character(err), fixed = TRUE),
  "external QC: any BUSCO-specific value requires genome mode"
)

merqury_only_qc <- data.frame(
  accession = "GCF_910000002.1",
  qc_busco_mode = "",
  merqury_qv = 48,
  merqury_completeness = 97,
  stringsAsFactors = FALSE
)
err <- try(
  rank_assembly_candidates(profile_meta, external_qc = merqury_only_qc),
  silent = TRUE
)
check(
  !inherits(err, "try-error"),
  "external QC: Merqury-only rows may omit BUSCO mode"
)

empty_evidence_qc <- data.frame(
  accession = "GCF_910000002.1",
  stringsAsFactors = FALSE
)
err <- try(
  rank_assembly_candidates(profile_meta, external_qc = empty_evidence_qc),
  silent = TRUE
)
check(
  inherits(err, "try-error") &&
    grepl("BUSCO or Merqury evidence", as.character(err), fixed = TRUE),
  "external QC: an accession-only row is rejected rather than treated as Merqury-only"
)

err <- try(
  rank_assembly_candidates(profile_meta, allow_qc_override = TRUE),
  silent = TRUE
)
check(
  inherits(err, "try-error") &&
    grepl("non-empty external QC", as.character(err), fixed = TRUE),
  "external QC: explicit override requires a QC input"
)

empty_qc <- data.frame(accession = character(), stringsAsFactors = FALSE)
err <- try(
  rank_assembly_candidates(
    profile_meta,
    external_qc = empty_qc,
    allow_qc_override = TRUE
  ),
  silent = TRUE
)
check(
  inherits(err, "try-error") &&
    grepl("non-empty external QC", as.character(err), fixed = TRUE),
  "external QC: explicit override rejects an empty QC table"
)

unversioned_qc <- strict_qc_row(
  "GCF_910000002", 96, 95, 1, 1, 3
)
err <- try(
  rank_assembly_candidates(profile_meta, external_qc = unversioned_qc),
  silent = TRUE
)
check(
  inherits(err, "try-error") &&
    grepl("exact versioned", as.character(err), fixed = TRUE),
  "external QC: unversioned accessions are rejected"
)

partly_unmatched_qc <- rbind(
  strict_qc_row("GCF_910000002.1", 96, 95, 1, 1, 3),
  strict_qc_row("GCF_910000002.2", 97, 96, 1, 1, 2),
  strict_qc_row("GCA_919999999.1", 97, 96, 1, 1, 2)
)
unmatched_warning <- ""
partly_matched <- withCallingHandlers(
  rank_assembly_candidates(profile_meta, external_qc = partly_unmatched_qc),
  warning = function(condition) {
    unmatched_warning <<- conditionMessage(condition)
    invokeRestart("muffleWarning")
  }
)
check(
  grepl("GCF_910000002.2", unmatched_warning, fixed = TRUE) &&
    grepl("GCA_919999999.1", unmatched_warning, fixed = TRUE) &&
    any(!is.na(partly_matched$audit$qc_busco_complete)),
  "external QC: version mismatches and typos are diagnosed while matched evidence remains attached"
)

all_unmatched_qc <- strict_qc_row(
  "GCA_919999998.1", 97, 96, 1, 1, 2
)
err <- try(suppressWarnings(rank_assembly_candidates(
  profile_meta,
  external_qc = all_unmatched_qc,
  allow_qc_override = TRUE
)), silent = TRUE)
check(
  inherits(err, "try-error") &&
    grepl("no accession matched", as.character(err), fixed = TRUE),
  "external QC: explicit override fails when no QC accession matches metadata"
)

derived_meta <- strict_pair(
  "Euk derived", "GCF_910000010.1", "GCA_910000011.1"
)
derived_qc <- rbind(
  strict_qc_row(
    "GCF_910000010.1", 96, duplicated = 2, fragmented = 1, missing = 3
  ),
  strict_qc_row(
    "GCA_910000011.1", 97, duplicated = 1, fragmented = 1, missing = 2
  )
)
derived <- rank_assembly_candidates(derived_meta, external_qc = derived_qc)$audit
check(identical(derived$qc_busco_single, c(94, 96)),
      "external QC: Single-copy is derived as Complete minus Duplicated")

invalid_qc_cases <- list(
  list(
    label = "non-numeric percentage",
    qc = within(
      strict_qc_row("GCF_910000020.1", 95, 93, 2, 2, 3),
      qc_busco_single <- "not-a-number"
    ),
    field = "qc_busco_single"
  ),
  list(
    label = "percentage above 100",
    qc = strict_qc_row("GCF_910000020.1", 101, 99, 2, 0, 0),
    field = "qc_busco_complete"
  ),
  list(
    label = "percentage below 0",
    qc = strict_qc_row("GCF_910000020.1", 99, 100, -1, 0, 1),
    field = "qc_busco_duplicated"
  ),
  list(
    label = "Complete below Duplicated",
    qc = strict_qc_row("GCF_910000020.1", 40, 0, 50, 10, 50),
    field = "qc_busco_complete"
  ),
  list(
    label = "Complete-Single-Duplicated inconsistency",
    qc = strict_qc_row("GCF_910000020.1", 95, 90, 2, 2, 3),
    field = "qc_busco_single"
  ),
  list(
    label = "Complete-Fragmented-Missing inconsistency",
    qc = strict_qc_row("GCF_910000020.1", 95, 93, 2, 1, 1),
    field = "qc_busco_complete"
  ),
  list(
    label = "internal-stop percentage above 100",
    qc = strict_qc_row("GCF_910000020.1", 95, 93, 2, 2, 3, 101),
    field = "qc_busco_internal_stop"
  )
)
validation_meta <- row("GCF_910000020.1", "Euk validation", reference = TRUE,
                       ncbi_busco_complete = 0.95)
for (case in invalid_qc_cases) {
  err <- try(rank_assembly_candidates(validation_meta, external_qc = case$qc),
             silent = TRUE)
  check(inherits(err, "try-error") &&
          grepl("GCF_910000020.1", as.character(err), fixed = TRUE) &&
          grepl(case$field, as.character(err), fixed = TRUE),
        paste("external QC validation:", case$label,
              "identifies accession and field"))
}

rounding_boundary_qc <- strict_qc_row(
  "GCF_910000020.1", 95, 92.8, 2, 2.2, 3
)
err <- try(rank_assembly_candidates(
  validation_meta, external_qc = rounding_boundary_qc
), silent = TRUE)
check(!inherits(err, "try-error"),
      "external QC validation: an exact 0.2 BUSCO rounding difference is accepted")

legacy_meta <- strict_pair(
  "Euk legacy", "GCF_910000030.1", "GCA_910000031.1"
)
legacy_qc <- rbind(
  strict_qc_row("GCF_910000030.1", 95, fragmented = 2, missing = 3),
  strict_qc_row("GCA_910000031.1", 98, fragmented = 1, missing = 1)
)
legacy <- rank_assembly_candidates(
  legacy_meta, external_qc = legacy_qc, allow_qc_override = TRUE
)
legacy_reason <- paste(legacy$audit$review_reason, collapse = ";")
check(identical(legacy$best$accession, "GCF_910000030.1") &&
        grepl("incomplete_busco_comparison", legacy_reason, fixed = TRUE),
      "override: incomplete legacy BUSCO rows are accepted for audit but cannot override")

plectropomus_meta <- strict_pair(
  "Plectropomus leopardus", "GCF_910000040.1", "GCA_910000041.1"
)
plectropomus_qc <- rbind(
  strict_qc_row("GCF_910000040.1", 96, 94, 2, 1, 3),
  strict_qc_row("GCA_910000041.1", 98, 97, 1, 0.5, 1.5)
)
plectropomus_default <- rank_assembly_candidates(
  plectropomus_meta, external_qc = plectropomus_qc
)
plectropomus_override <- rank_assembly_candidates(
  plectropomus_meta, external_qc = plectropomus_qc,
  allow_qc_override = TRUE
)
plectropomus_audit <- plectropomus_override$audit
check(identical(plectropomus_default$best$accession, "GCF_910000040.1") &&
        any(plectropomus_default$audit$qc_preferred[
          plectropomus_default$audit$accession == "GCA_910000041.1"
        ]),
      "review: a distinct higher-Single-copy candidate is QC-preferred without changing the default")
check(identical(plectropomus_override$best$accession, "GCA_910000041.1") &&
        isTRUE(plectropomus_audit$override_applied[
          plectropomus_audit$accession == "GCA_910000041.1"
        ]) &&
        identical(plectropomus_override$best$selection_basis,
                  "explicit_busco_override"),
      "override: explicit flag permits strict BUSCO dominance by a distinct candidate")
check(
  identical(
    plectropomus_audit$baseline_rank[
      plectropomus_audit$accession == "GCF_910000040.1"
    ],
    1L
  ) &&
    identical(
      plectropomus_audit$final_rank[
        plectropomus_audit$accession == "GCA_910000041.1"
      ],
      1L
    ) &&
    identical(plectropomus_override$best$final_rank, 1L),
  "override: baseline rank is preserved while the selected override receives final rank 1"
)
check(any(grepl(
  "alternative_higher_single_copy",
  plectropomus_default$audit$review_reason,
  fixed = TRUE
)), "review: strict higher-Single-copy evidence records a stable reason code")

larimichthys_meta <- strict_pair(
  "Larimichthys crocea", "GCF_910000050.1", "GCA_910000051.1"
)
larimichthys_qc <- rbind(
  strict_qc_row("GCF_910000050.1", 95, 94, 1, 2, 3),
  strict_qc_row("GCA_910000051.1", 97, 95, 2, 1, 2)
)
larimichthys <- rank_assembly_candidates(
  larimichthys_meta, external_qc = larimichthys_qc,
  allow_qc_override = TRUE
)
check(identical(larimichthys$best$accession, "GCF_910000050.1") && any(grepl(
  "complete_gain_duplication_confounded",
  larimichthys$audit$review_reason,
  fixed = TRUE
)) && identical(
  larimichthys$audit$override_block_reason[
    larimichthys$audit$accession == "GCA_910000051.1"
  ],
  "higher_duplicated"
), "override: a Complete gain with worse duplication records its stable blocker")

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
check(
  isTRUE(takifugu$audit$qc_preferred[
    takifugu$audit$accession == "GCF_000180615.1"
  ]) &&
    identical(takifugu$audit$qc_preferred[
      takifugu$audit$accession == "GCA_000180615.2"
    ], FALSE),
  "paired accessions: QC preference is assigned to the baseline representative, not its alias"
)

one_sided_meta <- strict_pair(
  "Euk one-sided", "GCF_910000060.1", "GCA_910000061.1"
)
one_sided_qc <- strict_qc_row(
  "GCA_910000061.1", 98, 97, 1, 0.5, 1.5
)
one_sided <- rank_assembly_candidates(
  one_sided_meta, external_qc = one_sided_qc, allow_qc_override = TRUE
)
check(identical(one_sided$best$accession, "GCF_910000060.1") && any(grepl(
  "one_sided_external_qc", one_sided$audit$review_reason, fixed = TRUE
)), "review: one-sided external QC is recorded and cannot override")

merqury_one_sided_meta <- strict_pair(
  "Euk Merqury one-sided", "GCF_910000062.1", "GCA_910000063.1"
)
merqury_one_sided_qc <- data.frame(
  accession = "GCA_910000063.1",
  qc_busco_mode = "",
  merqury_qv = 52,
  merqury_completeness = 99,
  stringsAsFactors = FALSE
)
merqury_one_sided <- rank_assembly_candidates(
  merqury_one_sided_meta,
  external_qc = merqury_one_sided_qc,
  allow_qc_override = TRUE
)
check(
  identical(merqury_one_sided$best$accession, "GCF_910000062.1") &&
    any(grepl(
      "one_sided_external_qc",
      merqury_one_sided$audit$review_reason,
      fixed = TRUE
    )) &&
    identical(
      sort(unique(merqury_one_sided$review$accession)),
      c("GCA_910000063.1", "GCF_910000062.1")
    ),
  "review: one-sided Merqury-only evidence reviews every candidate without overriding"
)

warning_meta <- rbind(
  strict_pair(
    "Euk low-complete", "GCF_910000070.1", "GCA_910000071.1"
  ),
  strict_pair(
    "Euk high-duplicate", "GCF_910000072.1", "GCA_910000073.1"
  )
)
warning_qc <- do.call(rbind, list(
  strict_qc_row("GCF_910000070.1", 89, 87, 2, 3, 8),
  strict_qc_row("GCA_910000071.1", 88, 86, 2, 3, 9),
  strict_qc_row("GCF_910000072.1", 96, 90, 6, 1, 3),
  strict_qc_row("GCA_910000073.1", 95, 90, 5, 1, 4)
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
  "GCF_910000074.1", "Euk solo warning", reference = TRUE,
  ncbi_busco_complete = 0.89
)
single_warning_qc <- strict_qc_row("GCF_910000074.1", 89, 87, 2, 3, 8)
single_warning_audit <- rank_assembly_candidates(
  single_warning_meta, external_qc = single_warning_qc
)$audit
check(grepl(
  "baseline_busco_complete_below_90",
  single_warning_audit$review_reason,
  fixed = TRUE
), "review: baseline BUSCO warnings do not require a distinct alternative")

qc_order_meta <- do.call(rbind, list(
  strict_pair("Euk QC D", "GCF_910000080.1", "GCA_910000081.1"),
  strict_pair("Euk QC M", "GCF_910000082.1", "GCA_910000083.1"),
  strict_pair("Euk QC F", "GCF_910000084.1", "GCA_910000085.1"),
  strict_pair("Euk QC metadata", "GCF_910000086.1", "GCA_910000087.1")
))
qc_order_qc <- do.call(rbind, list(
  strict_qc_row("GCF_910000080.1", 95, 90, 5, 2, 3),
  strict_qc_row("GCA_910000081.1", 94, 90, 4, 3, 3),
  strict_qc_row("GCF_910000082.1", 94, 90, 4, 2, 4),
  strict_qc_row("GCA_910000083.1", 94, 90, 4, 3, 3),
  strict_qc_row("GCF_910000084.1", 94, 90, 4, 2, 4),
  strict_qc_row("GCA_910000085.1", 94, 90, 4, 1.9, 4),
  strict_qc_row("GCF_910000086.1", 95, 92, 3, 2, 3),
  strict_qc_row("GCA_910000087.1", 95, 92, 3, 2, 3)
))
qc_order_audit <- rank_assembly_candidates(
  qc_order_meta, shortlist_k = Inf, external_qc = qc_order_qc
)$audit
qc_pick <- function(species) {
  qc_order_audit$accession[
    qc_order_audit$species == species & qc_order_audit$qc_preferred
  ]
}
check(identical(qc_pick("Euk QC D"), "GCA_910000081.1"),
      "QC preference: lower Duplicated breaks a Single-copy tie")
check(identical(qc_pick("Euk QC M"), "GCA_910000083.1"),
      "QC preference: lower Missing follows Duplicated")
check(identical(qc_pick("Euk QC F"), "GCA_910000085.1"),
      "QC preference: lower Fragmented follows Missing")
check(identical(qc_pick("Euk QC metadata"), "GCF_910000086.1"),
      "QC preference: existing metadata order is the final tie-breaker")

strict_blocker_meta <- rbind(
  strict_pair(
    "Euk fragmented blocker", "GCF_910000090.1", "GCA_910000091.1"
  ),
  strict_pair(
    "Euk missing blocker", "GCF_910000092.1", "GCA_910000093.1"
  )
)
strict_blocker_qc <- do.call(rbind, list(
  strict_qc_row("GCF_910000090.1", 95, 94, 1, 2, 3),
  strict_qc_row("GCA_910000091.1", 96, 95, 1, 3, 1),
  strict_qc_row("GCF_910000092.1", 95, 94, 1, 2, 3),
  strict_qc_row("GCA_910000093.1", 96, 95, 1, 0, 4)
))
strict_blockers <- rank_assembly_candidates(
  strict_blocker_meta, external_qc = strict_blocker_qc,
  allow_qc_override = TRUE
)
check(identical(
  strict_blockers$best$accession[
    strict_blockers$best$species == "Euk fragmented blocker"
  ],
  "GCF_910000090.1"
) && identical(
  strict_blockers$audit$override_block_reason[
    strict_blockers$audit$accession == "GCA_910000091.1"
  ],
  "higher_fragmented"
), "override: worse Fragmented records a stable blocker code")
check(identical(
  strict_blockers$best$accession[
    strict_blockers$best$species == "Euk missing blocker"
  ],
  "GCF_910000092.1"
) && identical(
  strict_blockers$audit$override_block_reason[
    strict_blockers$audit$accession == "GCA_910000093.1"
  ],
  "higher_missing"
), "override: worse Missing records a stable blocker code")

multi_override_meta <- rbind(
  strict_pair(
    "Euk multiple overrides", "GCF_910000100.1", "GCA_910000101.1"
  ),
  row(
    "GCA_910000102.1", "Euk multiple overrides", annotated = FALSE,
    ncbi_busco_complete = 0.97, contig_n50 = 4e6
  )
)
multi_override_qc <- rbind(
  strict_qc_row("GCF_910000100.1", 95, 94, 1, 2, 3),
  strict_qc_row("GCA_910000101.1", 96, 95, 1, 1.5, 2.5),
  strict_qc_row("GCA_910000102.1", 97, 96, 1, 1, 2)
)
multi_override <- rank_assembly_candidates(
  multi_override_meta, external_qc = multi_override_qc,
  allow_qc_override = TRUE
)
check(
  identical(multi_override$best$accession, "GCA_910000102.1") &&
    identical(
      multi_override$audit$override_block_reason[
        multi_override$audit$accession == "GCA_910000101.1"
      ],
      "not_top_override_candidate"
    ),
  "override: the strongest valid candidate wins and other valid candidates are recorded"
)

multi_blocker_meta <- strict_pair(
  "Euk multiple blockers", "GCF_910000110.1", "GCA_910000111.1"
)
multi_blocker_qc <- rbind(
  strict_qc_row("GCF_910000110.1", 95, 94, 1, 2, 3),
  strict_qc_row("GCA_910000111.1", 96, 93, 3, 3, 1)
)
multi_blocker <- rank_assembly_candidates(
  multi_blocker_meta, external_qc = multi_blocker_qc,
  allow_qc_override = TRUE
)
check(
  identical(
    multi_blocker$audit$override_block_reason[
      multi_blocker$audit$accession == "GCA_910000111.1"
    ],
    "not_higher_single_copy;higher_duplicated;higher_fragmented"
  ),
  "override: multiple blocker codes retain their stable evaluation order"
)

if (fail > 0L) { cat("\n", fail, " test(s) FAILED\n", sep = ""); quit(status = 1) }
cat("\nAll Stage 0 quality-selection tests passed.\n")
