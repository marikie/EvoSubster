#!/usr/bin/env Rscript

# Behavior tests for the two-phase Stage 0 assembly selection:
#   1. hard quality exclusions and reference-anchored shortlisting;
#   2. prokaryote CheckM/ANI ranking or eukaryote external BUSCO-genome ranking.
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

check(identical(pick("Prok alpha"), "P_best"),
      "prokaryote: higher CheckM completeness and lower contamination can beat reference")
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
check(is.data.frame(ranked$review) && nrow(ranked$review) == 0L,
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

if (fail > 0L) { cat("\n", fail, " test(s) FAILED\n", sep = ""); quit(status = 1) }
cat("\nAll Stage 0 quality-selection tests passed.\n")
