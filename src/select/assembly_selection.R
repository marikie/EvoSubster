# Pure Stage-0 assembly screening and ranking helpers.
# This module intentionally uses base R only so quality rules can be tested
# without loading phylogeny packages or contacting NCBI.

STAGE0_ASSEMBLY_LEVEL_RANK <- c(
  "Complete Genome" = 1,
  "Chromosome" = 2,
  "Scaffold" = 3,
  "Contig" = 4
)

STAGE0_METADATA_COLUMNS <- c(
  "source_database", "paired_accession", "assembly_status", "is_atypical",
  "atypical_warnings", "ani_check_status",
  "assembly_type", "refseq_category", "assembly_level", "contig_n50",
  "total_sequence_length", "total_ungapped_length", "number_of_contigs",
  "number_of_scaffolds", "scaffold_n50", "checkm_completeness",
  "checkm_contamination", "busco_lineage", "busco_version", "busco_complete",
  "busco_fragmented", "busco_missing", "has_annotation", "release_date"
)

STAGE0_EXTERNAL_QC_COLUMNS <- c(
  "qc_busco_mode", "qc_busco_lineage", "qc_busco_version",
  "qc_busco_complete", "qc_busco_single", "qc_busco_duplicated",
  "qc_busco_fragmented", "qc_busco_missing", "qc_busco_internal_stop",
  "merqury_qv", "merqury_completeness"
)

STAGE0_BUSCO_COMPONENT_COLUMNS <- c(
  "qc_busco_complete", "qc_busco_single", "qc_busco_duplicated",
  "qc_busco_fragmented", "qc_busco_missing"
)

STAGE0_BUSCO_ROUNDING_TOLERANCE <- 0.2 + sqrt(.Machine$double.eps)

stage0_add_missing_columns <- function(df, columns, value = NA) {
  for (column in setdiff(columns, names(df))) df[[column]] <- value
  df
}

stage0_numeric <- function(x) suppressWarnings(as.numeric(x))

stage0_true <- function(x) tolower(trimws(as.character(x))) %in% c("true", "1", "yes")

stage0_nonempty <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))

stage0_level_rank <- function(level) {
  rank <- unname(STAGE0_ASSEMBLY_LEVEL_RANK[as.character(level)])
  rank[is.na(rank)] <- length(STAGE0_ASSEMBLY_LEVEL_RANK) + 1L
  as.integer(rank)
}

stage0_primary_type_rank <- function(assembly_type) {
  type <- tolower(trimws(as.character(assembly_type)))
  rank <- rep(2L, length(type))
  rank[type %in% c("haploid", "primary", "primary assembly")] <- 1L
  rank[!nzchar(type) | is.na(type)] <- 3L
  rank
}

stage0_attach_external_qc <- function(meta, external_qc = NULL) {
  meta <- stage0_add_missing_columns(meta, STAGE0_EXTERNAL_QC_COLUMNS)
  if (is.null(external_qc) || !nrow(external_qc)) return(meta)
  if (!"accession" %in% names(external_qc)) {
    stop("External assembly QC must contain an accession column.", call. = FALSE)
  }
  if (anyDuplicated(external_qc$accession)) {
    stop("External assembly QC contains duplicate accessions.", call. = FALSE)
  }

  external_qc <- stage0_add_missing_columns(external_qc, STAGE0_EXTERNAL_QC_COLUMNS)
  invalid_mode <- stage0_nonempty(external_qc$qc_busco_mode) &
    tolower(external_qc$qc_busco_mode) != "genome"
  if (any(invalid_mode)) {
    stop("External BUSCO results must use qc_busco_mode=genome.", call. = FALSE)
  }

  match_index <- match(meta$accession, external_qc$accession)
  for (column in STAGE0_EXTERNAL_QC_COLUMNS) {
    meta[[column]] <- external_qc[[column]][match_index]
  }
  meta
}

stage0_validate_comparable_busco <- function(meta) {
  species_with_qc <- unique(meta$species[
    tolower(meta$qc_busco_mode) == "genome" & !is.na(meta$qc_busco_mode)
  ])
  for (species in species_with_qc) {
    rows <- meta$species == species & tolower(meta$qc_busco_mode) == "genome"
    rows[is.na(rows)] <- FALSE
    if (any(!stage0_nonempty(meta$qc_busco_lineage[rows])) ||
        any(!stage0_nonempty(meta$qc_busco_version[rows]))) {
      stop("External QC candidates for ", species,
           " must record BUSCO lineage and version.", call. = FALSE)
    }
    lineages <- unique(meta$qc_busco_lineage[rows & stage0_nonempty(meta$qc_busco_lineage)])
    versions <- unique(meta$qc_busco_version[rows & stage0_nonempty(meta$qc_busco_version)])
    if (length(lineages) > 1L || length(versions) > 1L) {
      stop("External QC candidates for ", species,
           " must use the same BUSCO lineage and version.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

stage0_busco_error <- function(meta, row, field, detail) {
  stop(
    "Invalid external BUSCO value for accession ", meta$accession[row],
    " field ", field, ": ", detail, ".",
    call. = FALSE
  )
}

stage0_validate_and_derive_busco <- function(meta) {
  percentage_columns <- c(
    STAGE0_BUSCO_COMPONENT_COLUMNS, "qc_busco_internal_stop"
  )
  for (field in percentage_columns) {
    invalid <- !is.na(meta[[field]]) &
      (!is.finite(meta[[field]]) | meta[[field]] < 0 | meta[[field]] > 100)
    if (any(invalid)) {
      stage0_busco_error(
        meta, which(invalid)[1], field, "percentage must be in [0, 100]"
      )
    }
  }

  complete_below_duplicated <- !is.na(meta$qc_busco_complete) &
    !is.na(meta$qc_busco_duplicated) &
    meta$qc_busco_complete < meta$qc_busco_duplicated
  if (any(complete_below_duplicated)) {
    stage0_busco_error(
      meta, which(complete_below_duplicated)[1], "qc_busco_complete",
      "Complete must be at least Duplicated"
    )
  }

  single_inconsistent <- !is.na(meta$qc_busco_complete) &
    !is.na(meta$qc_busco_single) & !is.na(meta$qc_busco_duplicated) &
    abs(meta$qc_busco_complete - meta$qc_busco_single -
          meta$qc_busco_duplicated) > STAGE0_BUSCO_ROUNDING_TOLERANCE
  if (any(single_inconsistent)) {
    stage0_busco_error(
      meta, which(single_inconsistent)[1], "qc_busco_single",
      "Complete must equal Single-copy plus Duplicated within 0.2"
    )
  }

  total_inconsistent <- !is.na(meta$qc_busco_complete) &
    !is.na(meta$qc_busco_fragmented) & !is.na(meta$qc_busco_missing) &
    abs(meta$qc_busco_complete + meta$qc_busco_fragmented +
          meta$qc_busco_missing - 100) > STAGE0_BUSCO_ROUNDING_TOLERANCE
  if (any(total_inconsistent)) {
    stage0_busco_error(
      meta, which(total_inconsistent)[1], "qc_busco_complete",
      "Complete plus Fragmented plus Missing must equal 100 within 0.2"
    )
  }

  derive_single <- is.na(meta$qc_busco_single) &
    !is.na(meta$qc_busco_complete) & !is.na(meta$qc_busco_duplicated)
  meta$qc_busco_single[derive_single] <-
    meta$qc_busco_complete[derive_single] - meta$qc_busco_duplicated[derive_single]
  meta
}

stage0_complete_busco <- function(df) {
  components_complete <- Reduce(
    `&`, lapply(STAGE0_BUSCO_COMPONENT_COLUMNS, function(field) {
      !is.na(df[[field]])
    })
  )
  tolower(trimws(as.character(df$qc_busco_mode))) == "genome" &
    stage0_nonempty(df$qc_busco_lineage) &
    stage0_nonempty(df$qc_busco_version) & components_complete
}

stage0_same_busco_run <- function(df, left, right) {
  all(tolower(trimws(as.character(df$qc_busco_mode[c(left, right)]))) ==
        "genome") &&
    identical(
      as.character(df$qc_busco_lineage[left]),
      as.character(df$qc_busco_lineage[right])
    ) &&
    identical(
      as.character(df$qc_busco_version[left]),
      as.character(df$qc_busco_version[right])
    )
}

stage0_qc_order <- function(df) {
  metadata_order <- stage0_baseline_order(df)$order
  metadata_rank <- integer(nrow(df))
  metadata_rank[metadata_order] <- seq_along(metadata_order)
  order(
    -df$qc_busco_single, df$qc_busco_duplicated, df$qc_busco_missing,
    df$qc_busco_fragmented, metadata_rank, na.last = TRUE
  )
}

stage0_add_review_reason <- function(candidates, row, reason) {
  existing <- candidates$review_reason[row]
  reasons <- if (nzchar(existing)) strsplit(existing, ";", fixed = TRUE)[[1]] else
    character()
  if (!reason %in% reasons) {
    candidates$review_reason[row] <- paste(c(reasons, reason), collapse = ";")
  }
  candidates$review_required[row] <- TRUE
  candidates
}

stage0_prepare_candidates <- function(meta, min_rel_contig_n50 = 0, external_qc = NULL) {
  if (!all(c("accession", "species") %in% names(meta))) {
    stop("Assembly metadata must contain accession and species columns.", call. = FALSE)
  }
  if (!is.finite(min_rel_contig_n50) || min_rel_contig_n50 < 0 ||
      min_rel_contig_n50 > 1) {
    stop("min_rel_contig_n50 must be in [0, 1].", call. = FALSE)
  }

  meta <- stage0_add_missing_columns(meta, STAGE0_METADATA_COLUMNS)
  meta <- stage0_attach_external_qc(meta, external_qc)
  stage0_validate_comparable_busco(meta)

  numeric_columns <- c(
    "contig_n50", "total_sequence_length", "total_ungapped_length",
    "number_of_contigs", "number_of_scaffolds", "scaffold_n50",
    "checkm_completeness", "checkm_contamination", "busco_complete",
    "busco_fragmented", "busco_missing", STAGE0_BUSCO_COMPONENT_COLUMNS,
    "qc_busco_internal_stop", "merqury_qv", "merqury_completeness"
  )
  for (column in numeric_columns) {
    raw_value <- meta[[column]]
    numeric_value <- stage0_numeric(raw_value)
    if (column %in% c(STAGE0_BUSCO_COMPONENT_COLUMNS,
                      "qc_busco_internal_stop")) {
      invalid_numeric <- stage0_nonempty(raw_value) & is.na(numeric_value)
      if (any(invalid_numeric)) {
        stage0_busco_error(
          meta, which(invalid_numeric)[1], column, "value must be numeric"
        )
      }
    }
    meta[[column]] <- numeric_value
  }
  meta <- stage0_validate_and_derive_busco(meta)

  meta$rel_contig_n50 <- ifelse(
    is.na(meta$total_ungapped_length) | meta$total_ungapped_length <= 0 |
      is.na(meta$contig_n50),
    0,
    meta$contig_n50 / meta$total_ungapped_length
  )
  meta$gap_fraction <- ifelse(
    is.na(meta$total_sequence_length) | meta$total_sequence_length <= 0 |
      is.na(meta$total_ungapped_length),
    NA_real_,
    pmax(0, 1 - meta$total_ungapped_length / meta$total_sequence_length)
  )
  meta$is_reference <- tolower(trimws(as.character(meta$refseq_category))) ==
    "reference genome"
  source_is_refseq <- toupper(trimws(as.character(meta$source_database))) ==
    "SOURCE_DATABASE_REFSEQ"
  source_is_refseq[is.na(source_is_refseq)] <- FALSE
  meta$is_refseq <- source_is_refseq | grepl("^GCF_", meta$accession)
  meta$assembly_level_rank <- stage0_level_rank(meta$assembly_level)
  meta$primary_type_rank <- stage0_primary_type_rank(meta$assembly_type)

  status_current <- tolower(trimws(as.character(meta$assembly_status))) == "current"
  atypical <- stage0_true(meta$is_atypical) | stage0_nonempty(meta$atypical_warnings)
  ani_failed <- toupper(trimws(as.character(meta$ani_check_status))) == "FAILED"
  type_text <- tolower(trimws(as.character(meta$assembly_type)))
  unsupported_type <- grepl(
    "hybrid|unresolved[ _-]*diploid|alternate[ _-]*haplotype",
    type_text
  )
  below_gate <- min_rel_contig_n50 > 0 & meta$rel_contig_n50 < min_rel_contig_n50

  meta$eligible <- TRUE
  meta$exclusion_reason <- ""
  set_reason <- function(mask, reason) {
    mask[is.na(mask)] <- FALSE
    targets <- meta$eligible & mask
    meta$eligible[targets] <<- FALSE
    meta$exclusion_reason[targets] <<- reason[targets]
  }
  set_reason(!status_current, rep("not_current", nrow(meta)))
  atypical_reason <- ifelse(
    stage0_nonempty(meta$atypical_warnings),
    paste0("atypical: ", meta$atypical_warnings),
    "atypical"
  )
  set_reason(atypical, atypical_reason)
  set_reason(ani_failed, rep("ani_failed", nrow(meta)))
  set_reason(
    unsupported_type,
    paste0("unsupported_assembly_type: ", meta$assembly_type)
  )
  set_reason(below_gate, rep("below_relative_contig_n50_gate", nrow(meta)))

  meta$selection_profile <- "general"
  for (species in unique(meta$species)) {
    rows <- meta$species == species
    has_prokaryote_qc <- any(!is.na(meta$checkm_completeness[rows])) ||
      any(stage0_nonempty(meta$ani_check_status[rows]))
    has_eukaryote_qc <- any(!is.na(meta$busco_complete[rows])) ||
      any(stage0_nonempty(meta$busco_lineage[rows])) ||
      any(tolower(meta$qc_busco_mode[rows]) == "genome", na.rm = TRUE)
    if (has_prokaryote_qc) meta$selection_profile[rows] <- "prokaryote"
    else if (has_eukaryote_qc) meta$selection_profile[rows] <- "eukaryote"
  }

  meta
}

stage0_assembly_equivalence_key <- function(accession, paired_accession) {
  vapply(seq_along(accession), function(i) {
    accession_value <- trimws(as.character(accession[i]))
    paired_value <- trimws(as.character(paired_accession[i]))
    if (!stage0_nonempty(paired_value)) return(accession_value)
    paste(sort(c(accession_value, paired_value)), collapse = "|")
  }, character(1))
}

stage0_baseline_order <- function(df) {
  if (df$selection_profile[1] == "eukaryote") {
    annotated <- stage0_true(df$has_annotation)
    tier <- rep(5L, nrow(df))
    tier[annotated] <- 4L
    tier[annotated & df$is_refseq] <- 3L
    tier[df$is_reference] <- 2L
    tier[df$is_reference & annotated] <- 1L
    return(list(
      order = order(
        tier, df$assembly_level_rank, -df$contig_n50, df$gap_fraction,
        df$primary_type_rank, df$accession, na.last = TRUE
      ),
      basis = "reference_metadata_baseline"
    ))
  }

  stage0_final_order(df)
}

stage0_phase1_order <- function(df) {
  profile <- df$selection_profile[1]
  ani_ok <- toupper(trimws(as.character(df$ani_check_status))) == "OK"
  common <- list(
    -as.integer(df$is_reference),
    df$assembly_level_rank,
    -df$contig_n50,
    df$number_of_contigs,
    df$gap_fraction,
    df$accession
  )
  if (profile == "prokaryote") {
    return(order(
      -as.integer(df$is_reference), -as.integer(ani_ok),
      -df$checkm_completeness, df$checkm_contamination,
      df$assembly_level_rank, df$number_of_contigs, -df$contig_n50,
      df$accession, na.last = TRUE
    ))
  }
  if (profile == "eukaryote") {
    return(stage0_baseline_order(df)$order)
  }
  do.call(order, c(common, list(na.last = TRUE)))
}

stage0_final_order <- function(df) {
  profile <- df$selection_profile[1]
  ani_ok <- toupper(trimws(as.character(df$ani_check_status))) == "OK"
  qc_available <- tolower(df$qc_busco_mode) == "genome" &
    !is.na(df$qc_busco_complete)

  if (profile == "prokaryote") {
    has_merqury <- any(!is.na(df$merqury_qv) | !is.na(df$merqury_completeness))
    return(list(
      order = order(
        -as.integer(ani_ok), -df$checkm_completeness, df$checkm_contamination,
        -df$merqury_completeness, -df$merqury_qv,
        df$assembly_level_rank, df$number_of_contigs, -df$contig_n50,
        df$gap_fraction, -as.integer(df$is_reference), df$accession,
        na.last = TRUE
      ),
      basis = if (has_merqury) "ani_checkm_merqury_assembly_quality" else
        "ani_checkm_assembly_quality"
    ))
  }

  if (profile == "eukaryote" && any(qc_available)) {
    return(list(
      order = order(
        -as.integer(qc_available), -df$qc_busco_complete,
        df$qc_busco_missing, df$qc_busco_fragmented,
        -df$merqury_completeness, -df$merqury_qv,
        df$assembly_level_rank, -df$contig_n50, df$gap_fraction,
        df$primary_type_rank, -as.integer(df$is_reference), df$accession,
        na.last = TRUE
      ),
      basis = "external_busco_genome_merqury"
    ))
  }

  if (profile == "eukaryote") {
    return(list(
      order = order(
        -as.integer(df$is_reference), df$assembly_level_rank,
        -df$contig_n50, df$gap_fraction, df$primary_type_rank,
        -as.integer(stage0_true(df$has_annotation)), df$accession,
        na.last = TRUE
      ),
      basis = "reference_metadata_provisional"
    ))
  }

  list(
    order = order(
      -as.integer(df$is_reference), df$assembly_level_rank,
      -df$contig_n50, df$number_of_contigs, df$gap_fraction,
      df$primary_type_rank, -as.integer(stage0_true(df$has_annotation)),
      df$accession, na.last = TRUE
    ),
    basis = "reference_metadata_general"
  )
}

rank_assembly_candidates <- function(meta, min_rel_contig_n50 = 0,
                                     shortlist_k = 3, external_qc = NULL,
                                     allow_qc_override = FALSE) {
  if (length(shortlist_k) != 1L || is.na(shortlist_k) || shortlist_k < 1) {
    stop("shortlist_k must be a positive integer or Inf.", call. = FALSE)
  }
  candidates <- stage0_prepare_candidates(meta, min_rel_contig_n50, external_qc)
  candidates$shortlist_rank <- NA_integer_
  candidates$shortlisted <- FALSE
  candidates$final_rank <- NA_integer_
  candidates$selected <- FALSE
  candidates$selection_basis <- ""
  candidates$baseline_selected <- FALSE
  candidates$qc_preferred <- FALSE
  candidates$review_required <- FALSE
  candidates$review_reason <- ""
  candidates$override_applied <- FALSE
  candidates$assembly_equivalence_key <- stage0_assembly_equivalence_key(
    candidates$accession, candidates$paired_accession
  )

  for (species in unique(candidates$species)) {
    eligible_index <- which(candidates$species == species & candidates$eligible)
    if (!length(eligible_index)) next

    eligible <- candidates[eligible_index, , drop = FALSE]
    phase1 <- stage0_phase1_order(eligible)
    ordered_index <- eligible_index[phase1]
    candidates$shortlist_rank[ordered_index] <- seq_along(ordered_index)

    keep_n <- if (is.infinite(shortlist_k)) length(ordered_index) else
      min(length(ordered_index), as.integer(shortlist_k))
    shortlist_index <- ordered_index[seq_len(keep_n)]
    candidates$shortlisted[shortlist_index] <- TRUE

    baseline <- stage0_baseline_order(
      candidates[shortlist_index, , drop = FALSE]
    )
    final_index <- shortlist_index[baseline$order]
    candidates$final_rank[final_index] <- seq_along(final_index)
    candidates$selection_basis[shortlist_index] <- baseline$basis
    candidates$baseline_selected[final_index[1]] <- TRUE
    candidates$selected[final_index[1]] <- TRUE

    if (eligible$selection_profile[1] == "eukaryote") {
      baseline_index <- final_index[1]
      complete_index <- shortlist_index[
        stage0_complete_busco(candidates[shortlist_index, , drop = FALSE])
      ]
      if (length(complete_index)) {
        qc_order <- stage0_qc_order(
          candidates[complete_index, , drop = FALSE]
        )
        candidates$qc_preferred[complete_index[qc_order[1]]] <- TRUE
      }

      distinct_alternatives <- setdiff(
        shortlist_index[
          candidates$assembly_equivalence_key[shortlist_index] !=
            candidates$assembly_equivalence_key[baseline_index]
        ],
        baseline_index
      )
      baseline_has_qc <- identical(
        tolower(trimws(as.character(candidates$qc_busco_mode[baseline_index]))),
        "genome"
      )
      baseline_complete <- baseline_index %in% complete_index

      if (baseline_has_qc &&
          !is.na(candidates$qc_busco_complete[baseline_index]) &&
          candidates$qc_busco_complete[baseline_index] < 90) {
        candidates <- stage0_add_review_reason(
          candidates, baseline_index, "baseline_busco_complete_below_90"
        )
      }
      if (baseline_has_qc &&
          !is.na(candidates$qc_busco_duplicated[baseline_index]) &&
          candidates$qc_busco_duplicated[baseline_index] > 5) {
        candidates <- stage0_add_review_reason(
          candidates, baseline_index, "baseline_busco_duplicated_above_5"
        )
      }

      if (length(distinct_alternatives)) {
        for (alternative_index in distinct_alternatives) {
          alternative_has_qc <- identical(
            tolower(trimws(as.character(
              candidates$qc_busco_mode[alternative_index]
            ))),
            "genome"
          )
          alternative_complete <- alternative_index %in% complete_index

          if (xor(baseline_has_qc, alternative_has_qc)) {
            candidates <- stage0_add_review_reason(
              candidates, alternative_index, "one_sided_external_qc"
            )
            candidates <- stage0_add_review_reason(
              candidates, alternative_index, "incomplete_busco_comparison"
            )
          } else if (baseline_has_qc && alternative_has_qc &&
                     (!baseline_complete || !alternative_complete)) {
            candidates <- stage0_add_review_reason(
              candidates, alternative_index, "incomplete_busco_comparison"
            )
          } else if (baseline_complete && alternative_complete &&
                     stage0_same_busco_run(
                       candidates, baseline_index, alternative_index
                     )) {
            if (candidates$qc_busco_single[alternative_index] >
                candidates$qc_busco_single[baseline_index]) {
              candidates <- stage0_add_review_reason(
                candidates, alternative_index,
                "alternative_higher_single_copy"
              )
            }
            if (candidates$qc_busco_complete[alternative_index] >
                  candidates$qc_busco_complete[baseline_index] &&
                candidates$qc_busco_duplicated[alternative_index] >
                  candidates$qc_busco_duplicated[baseline_index]) {
              candidates <- stage0_add_review_reason(
                candidates, alternative_index,
                "complete_gain_duplication_confounded"
              )
            }
          }
        }
      }

      if (isTRUE(allow_qc_override) && baseline_complete &&
          length(distinct_alternatives)) {
        override_candidates <- distinct_alternatives[vapply(
          distinct_alternatives,
          function(alternative_index) {
            alternative_index %in% complete_index &&
              stage0_same_busco_run(
                candidates, baseline_index, alternative_index
              ) &&
              candidates$qc_busco_single[alternative_index] >
                candidates$qc_busco_single[baseline_index] &&
              candidates$qc_busco_duplicated[alternative_index] <=
                candidates$qc_busco_duplicated[baseline_index] &&
              candidates$qc_busco_fragmented[alternative_index] <=
                candidates$qc_busco_fragmented[baseline_index] &&
              candidates$qc_busco_missing[alternative_index] <=
                candidates$qc_busco_missing[baseline_index]
          },
          logical(1)
        )]
        if (length(override_candidates)) {
          override_order <- stage0_qc_order(
            candidates[override_candidates, , drop = FALSE]
          )
          override_index <- override_candidates[override_order[1]]
          candidates$selected[baseline_index] <- FALSE
          candidates$selected[override_index] <- TRUE
          candidates$selection_basis[override_index] <-
            "explicit_busco_override"
          candidates$override_applied[override_index] <- TRUE
        }
      }
    }
  }

  shortlist <- candidates[candidates$eligible & candidates$shortlisted, , drop = FALSE]
  shortlist <- shortlist[order(shortlist$species, shortlist$shortlist_rank), , drop = FALSE]
  best <- candidates[candidates$eligible & candidates$selected, , drop = FALSE]
  best <- best[order(best$species), , drop = FALSE]
  review_species <- unique(candidates$species[candidates$review_required])
  review <- candidates[candidates$species %in% review_species, , drop = FALSE]
  review <- review[order(review$species, review$final_rank), , drop = FALSE]
  rownames(candidates) <- rownames(shortlist) <- rownames(best) <- rownames(review) <- NULL

  list(audit = candidates, shortlist = shortlist, best = best, review = review)
}

select_best_assemblies <- function(meta, min_rel_contig_n50 = 0,
                                   shortlist_k = 3, external_qc = NULL,
                                   allow_qc_override = FALSE) {
  rank_assembly_candidates(
    meta,
    min_rel_contig_n50 = min_rel_contig_n50,
    shortlist_k = shortlist_k,
    external_qc = external_qc,
    allow_qc_override = allow_qc_override
  )$best
}
