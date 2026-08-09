#!/usr/bin/env Rscript

# Select substitution-analysis trios from a Newick tree.
#
# Original script written by Ryota Ishii; rewritten to apply the thesis
# (section 2.3) rule -- measured sequence percent identity plus genus names --
# instead of tree distances.
#
# One outgroup per ingroup pair.  The tree is used only for TOPOLOGY: for an
# ingroup couple it gives the candidate outgroups ordered nearest -> farthest
# (sister of their MRCA, then sister of the grandparent, ...).  The tree is
# ultrametric, so its distances cannot rank outgroups -- the real decision needs
# last-train identity.  For each pair we walk the candidates nearest-first,
# download + last-train on demand, and STOP at the first outgroup that passes the
# thesis rule.  The strict thesis ordering (idt_12 < idt_23 AND idt_13 < idt_23)
# is itself the "not too close" boundary, and --idt-threshold is the "not too far"
# boundary, so nearest-first + first-pass yields the closest still-external
# outgroup while training as few pairs as possible.
#
# Stage 0  prune the tree to one assembly per species      [NCBI metadata, no download]
# Stage 1  per ingroup pair: search outgroups near->far,   [dwl_organism.sh + last_train.sh
#          stop at the first thesis-rule pass -> one trio    on demand; trio_filter.py to judge]
#
# Usage:
#   Rscript src/select/trio_selection.R --tree Primates_tree.nwk [options]
#
# Depends on: ape, dplyr (R); python3, LAST, curl, jq (external).

SCRIPT_DIR <- local({
  sourced_file <- tryCatch(sys.frame(1)$ofile, error = function(...) NULL)
  if (!is.null(sourced_file) && nzchar(sourced_file)) {
    return(dirname(normalizePath(sourced_file)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
})
SRC_DIR <- dirname(SCRIPT_DIR)

FETCH_METADATA <- file.path(SCRIPT_DIR, "fetch_assembly_metadata.py")
TRIO_FILTER <- file.path(SCRIPT_DIR, "trio_filter.py")
DWL_ORGANISM <- file.path(SRC_DIR, "dwl_organism.sh")
LAST_TRAIN <- file.path(SRC_DIR, "align", "last_train.sh")
ASSEMBLY_SELECTION <- file.path(SCRIPT_DIR, "assembly_selection.R")
source(ASSEMBLY_SELECTION)

IDENTITY_PREFIX <- "# substitution percent identity:"

# --- options ----------------------------------------------------------------

defaults <- list(
  tree = NULL,
  idt_threshold = 80,
  max_outgroup_tries = 5,
  ingroup_pairing = "matching",
  min_rel_contig_n50 = 0,
  stage0_top_k = 3L,
  assembly_qc = NULL,
  allow_qc_override = FALSE,
  genome_dir = "./genomes",
  out_dir = "./results/trio_selection",
  train_cache_dir = NULL,
  date = format(Sys.Date(), "%Y%m%d"),
  threads = 8,
  dry_run = FALSE,
  keep_unused_species_data = FALSE
)

usage <- function() {
  cat(
    "Usage: Rscript src/select/trio_selection.R --tree <file.nwk> [options]",
    "",
    "Options:",
    "  --tree FILE             Newick tree; leaf labels are accession-free complete taxon names",
    "                          with spaces encoded as underscores (for example,",
    "                          Chaunax_sp._Z400). Legacy accession-suffixed trees must first",
    "                          be converted with src/select/strip_newick_accessions.py.",
    "  --idt-threshold N       Every pairwise identity must exceed this (default: 80).",
    "  --max-outgroup-tries N  Give up on an ingroup pair after training this many outgroup",
    "                          candidates without a thesis-rule pass (default: 5). Caps the",
    "                          per-pair last-train cost.",
    "  --ingroup-pairing MODE  Which ingroup couples to consider (default: matching).",
    "                            matching  greedy closest-first disjoint matching -- each",
    "                                      species used at most once, ~n/2 independent",
    "                                      sister-pair couples (phylogenetic-pairs design).",
    "                            all       every tip pair (exhaustive; C(n,2) couples).",
    "  --min-rel-contig-n50 F  Drop a species unless some current assembly has relative contig",
    "                          N50 (contig_n50 / total_ungapped_length) >= F. Size-normalized,",
    "                          so one threshold works across lineages -- absolute N50 is",
    "                          meaningless across taxa (a perfect ~12 Mb yeast reference has a",
    "                          0.92 Mb contig N50). Default 0 = off (keep every species);",
    "                          e.g. 0.005 = 0.5%.",
    "  --stage0-top-k N        Keep N candidates per species in the metadata shortlist for QC audit and optional explicit strict override (default: 3).",
    "                          The eligible NCBI reference",
    "                          genome, when present, is always anchored in this shortlist.",
    "  --assembly-qc FILE      Optional TSV with accession plus comparable external QC:",
    "                          BUSCO genome-mode fields (qc_busco_*) and optional Merqury",
    "                          fields. NCBI annotation BUSCO is recorded but is not treated",
    "                          as assembly-level BUSCO genome-mode evidence.",
    "  --allow-qc-override    Explicitly allow a strict BUSCO-based override of the metadata",
    "                          baseline. Off by default; external QC is otherwise audit-only.",
    "  --genome-dir DIR        Genome storage (default: ./genomes).",
    "  --out-dir DIR           Output directory (default: ./results/trio_selection).",
    "  --train-cache-dir DIR   Directory for cached last-train files (default: <out-dir>/train_cache).",
    "  --date YYYYMMDD         Run date, used in cache filenames (default: today).",
    "  --threads N             Threads for lastdb/last-train (default: 8).",
    "  --dry-run               Enumerate ingroup pairs + candidate-outgroup order only;",
    "                          download nothing and run no last-train.",
    "  --keep-unused-species-data",
    "                          Keep genome and train artifacts created for candidates that are",
    "                          absent from every selected trio (default: remove them).",
    "",
    sep = "\n"
  )
}

load_runtime_dependencies <- function() {
  packages <- c("ape", "dplyr")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      ". Install with install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), ")).",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages({
    library(ape)
    library(dplyr)
  })
}

parse_args <- function(argv) {
  opts <- defaults
  i <- 1
  while (i <= length(argv)) {
    key <- argv[i]
    take <- function() {
      if (i + 1 > length(argv)) stop("Missing value for ", key, call. = FALSE)
      argv[i + 1]
    }
    switch(key,
      "--tree"               = { opts$tree <- take(); i <- i + 2 },
      "--idt-threshold"      = { opts$idt_threshold <- as.numeric(take()); i <- i + 2 },
      "--max-outgroup-tries" = { opts$max_outgroup_tries <- as.integer(take()); i <- i + 2 },
      "--ingroup-pairing"    = { opts$ingroup_pairing <- take(); i <- i + 2 },
      "--min-rel-contig-n50" = { opts$min_rel_contig_n50 <- as.numeric(take()); i <- i + 2 },
      "--stage0-top-k"       = { opts$stage0_top_k <- as.integer(take()); i <- i + 2 },
      "--assembly-qc"        = { opts$assembly_qc <- take(); i <- i + 2 },
      "--allow-qc-override"  = { opts$allow_qc_override <- TRUE; i <- i + 1 },
      "--genome-dir"         = { opts$genome_dir <- take(); i <- i + 2 },
      "--out-dir"            = { opts$out_dir <- take(); i <- i + 2 },
      "--train-cache-dir"    = { opts$train_cache_dir <- take(); i <- i + 2 },
      "--date"               = { opts$date <- take(); i <- i + 2 },
      "--threads"            = { opts$threads <- as.integer(take()); i <- i + 2 },
      "--dry-run"            = { opts$dry_run <- TRUE; i <- i + 1 },
      "--keep-unused-species-data" = { opts$keep_unused_species_data <- TRUE; i <- i + 1 },
      "-h"                   = { usage(); quit(status = 0) },
      "--help"               = { usage(); quit(status = 0) },
      stop("Unknown option: ", key, call. = FALSE)
    )
  }
  if (is.null(opts$tree)) {
    usage()
    stop("--tree is required.", call. = FALSE)
  }
  if (!opts$ingroup_pairing %in% c("matching", "all")) {
    stop("--ingroup-pairing must be 'matching' or 'all' (got: ", opts$ingroup_pairing, ")",
         call. = FALSE)
  }
  if (is.na(opts$min_rel_contig_n50) || opts$min_rel_contig_n50 < 0 ||
      opts$min_rel_contig_n50 > 1) {
    stop("--min-rel-contig-n50 must be a fraction in [0, 1] (got: ", opts$min_rel_contig_n50,
         "); it is contig_n50 / total_ungapped_length, e.g. 0.005 for 0.5%.", call. = FALSE)
  }
  if (is.na(opts$stage0_top_k) || opts$stage0_top_k < 1) {
    stop("--stage0-top-k must be a positive integer (got: ", opts$stage0_top_k, ")",
         call. = FALSE)
  }
  opts
}

# --- helpers ----------------------------------------------------------------

log_stage <- function(...) cat("[trio_selection]", ..., "\n")

run <- function(command, args, quiet_stderr = FALSE) {
  output <- suppressWarnings(system2(command, args, stdout = TRUE,
                                     stderr = if (quiet_stderr) FALSE else ""))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop("Command failed (", status, "): ", command, " ", paste(args, collapse = " "),
         call. = FALSE)
  }
  output
}

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

write_tsv_atomic <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit({
    if (path_exists_without_following(tmp)) unlink(tmp)
  }, add = TRUE)
  write_tsv(df, tmp)
  if (!file.rename(tmp, path)) {
    stop("Could not atomically replace TSV file: ", path, call. = FALSE)
  }
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
             colClasses = "character")
}

empty_owned_artifacts <- function() {
  data.frame(
    accession_1 = character(),
    accession_2 = character(),
    storage = character(),
    artifact_type = character(),
    relative_path = character(),
    stringsAsFactors = FALSE
  )
}

read_artifact_manifest <- function(path, requested_accession) {
  if (!file.exists(path)) {
    stop("dwl_organism.sh did not create its artifact manifest for ", requested_accession,
         call. = FALSE)
  }
  manifest <- read_tsv(path)
  expected <- c("accession", "artifact_type", "relative_path")
  if (!identical(names(manifest), expected)) {
    stop("Invalid downloader artifact manifest header: ", path, call. = FALSE)
  }
  if (!nrow(manifest)) return(empty_owned_artifacts())

  allowed_types <- c("directory_tree", "directory", "file", "symlink")
  invalid <- is.na(manifest$accession) | !nzchar(manifest$accession) |
    manifest$accession != requested_accession |
    is.na(manifest$artifact_type) | !manifest$artifact_type %in% allowed_types |
    is.na(manifest$relative_path) | !nzchar(manifest$relative_path)
  if (any(invalid)) {
    stop("Invalid row in downloader artifact manifest: ", path, call. = FALSE)
  }

  data.frame(
    accession_1 = manifest$accession,
    accession_2 = "",
    storage = "genome",
    artifact_type = manifest$artifact_type,
    relative_path = manifest$relative_path,
    stringsAsFactors = FALSE
  )
}

path_exists_without_following <- function(path) {
  link_target <- Sys.readlink(path)
  file.exists(path) || dir.exists(path) || (!is.na(link_target) && nzchar(link_target))
}

prepare_cleanup_target <- function(relative_path, root) {
  if (is.na(relative_path) || !nzchar(relative_path) ||
      startsWith(relative_path, "/") ||
      grepl("^[A-Za-z]:[/\\\\]", relative_path)) {
    stop("Unsafe cleanup path: ", relative_path, call. = FALSE)
  }
  components <- strsplit(relative_path, "/", fixed = TRUE)[[1]]
  if (any(!nzchar(components)) || any(components %in% c(".", ".."))) {
    stop("Unsafe cleanup path: ", relative_path, call. = FALSE)
  }

  normalized_root <- normalizePath(root, mustWork = TRUE)
  normalized_parent <- normalizePath(
    file.path(normalized_root, dirname(relative_path)), mustWork = TRUE
  )
  root_prefix <- paste0(normalized_root, .Platform$file.sep)
  if (normalized_parent != normalized_root && !startsWith(normalized_parent, root_prefix)) {
    stop("Cleanup path escapes configured root: ", relative_path, call. = FALSE)
  }

  target <- file.path(normalized_parent, basename(relative_path))
  if (identical(target, normalized_root)) {
    stop("Refusing to delete a cleanup root", call. = FALSE)
  }
  target
}

cleanup_owned_artifacts <- function(artifacts, selected, genome_root, train_cache_root,
                                    audit_path) {
  required_columns <- c(
    "accession_1", "accession_2", "storage", "artifact_type", "relative_path"
  )
  if (!all(required_columns %in% names(artifacts))) {
    stop("Owned artifact table is missing required columns", call. = FALSE)
  }

  selected_accessions <- character()
  if (nrow(selected)) {
    selected_columns <- c("out_acc", "in1_acc", "in2_acc")
    if (!all(selected_columns %in% names(selected))) {
      stop("Selected trio table is missing accession columns", call. = FALSE)
    }
    selected_accessions <- unique(unname(unlist(selected[selected_columns], use.names = FALSE)))
    selected_accessions <- selected_accessions[!is.na(selected_accessions) &
                                                 nzchar(selected_accessions)]
  }

  audit <- data.frame(
    accessions = character(), artifact_type = character(), path = character(),
    action = character(), reason = character(), stringsAsFactors = FALSE
  )
  dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
  if (!nrow(artifacts)) {
    write_tsv_atomic(audit, audit_path)
    return(invisible(audit))
  }

  accessions <- vapply(seq_len(nrow(artifacts)), function(i) {
    values <- c(artifacts$accession_1[i], artifacts$accession_2[i])
    paste(values[!is.na(values) & nzchar(values)], collapse = ",")
  }, character(1))
  pair_artifact <- !is.na(artifacts$accession_2) & nzchar(artifacts$accession_2)
  retained <- artifacts$accession_1 %in% selected_accessions &
    (artifacts$storage == "genome" | !pair_artifact |
       artifacts$accession_2 %in% selected_accessions)
  selected_reason <- ifelse(
    artifacts$storage == "genome" | !pair_artifact,
    "selected_accession", "selected_accession_pair"
  )
  unused_reason <- ifelse(
    artifacts$storage == "genome" | !pair_artifact,
    "unused_accession", "unused_accession_pair"
  )
  audit <- data.frame(
    accessions = accessions,
    artifact_type = artifacts$artifact_type,
    path = artifacts$relative_path,
    action = ifelse(retained, "retained", "pending"),
    reason = ifelse(retained, selected_reason, unused_reason),
    stringsAsFactors = FALSE
  )
  # Replace a previous run's report before any destructive preflight begins.
  write_tsv_atomic(audit, audit_path)

  active_index <- NA_integer_
  tryCatch({
    allowed_types <- c("directory_tree", "directory", "file", "symlink")
    invalid_train_type <- artifacts$storage == "train_cache" &
      !artifacts$artifact_type %in% c("file", "directory_tree")
    if (any(!artifacts$storage %in% c("genome", "train_cache")) ||
        any(!artifacts$artifact_type %in% allowed_types) || any(invalid_train_type)) {
      stop("Owned artifact table contains an invalid storage or artifact type", call. = FALSE)
    }

    roots <- c(genome = genome_root, train_cache = train_cache_root)
    targets <- character(nrow(artifacts))
    for (i in seq_len(nrow(artifacts))) {
      active_index <- i
      targets[i] <- prepare_cleanup_target(
        artifacts$relative_path[i], roots[[artifacts$storage[i]]]
      )
    }
    active_index <- NA_integer_
    audit$path <- targets
    write_tsv_atomic(audit, audit_path)

    delete_indices <- which(!retained)
    if (length(delete_indices)) {
      duplicate_targets <- targets[duplicated(targets) | duplicated(targets, fromLast = TRUE)]
      if (length(duplicate_targets)) {
        duplicate_targets <- unique(duplicate_targets)
        conflict <- any(vapply(duplicate_targets, function(target) {
          states <- retained[targets == target]
          any(states) && any(!states)
        }, logical(1)))
        if (conflict) {
          stop("Cleanup target is both retained and scheduled for deletion", call. = FALSE)
        }
        stop("Owned artifact table contains a duplicate cleanup target", call. = FALSE)
      }

      recursive_targets <- targets[
        !retained & artifacts$artifact_type == "directory_tree"
      ]
      retained_targets <- targets[retained]
      for (target in recursive_targets) {
        prefix <- paste0(target, .Platform$file.sep)
        if (any(retained_targets == target | startsWith(retained_targets, prefix))) {
          stop("Recursive cleanup target contains a retained artifact: ", target,
               call. = FALSE)
        }
      }

      # Validate every target before deleting any of them.
      for (i in delete_indices) {
        active_index <- i
        target <- targets[i]
        link_target <- Sys.readlink(target)
        is_link <- !is.na(link_target) && nzchar(link_target)
        valid_type <- switch(
          artifacts$artifact_type[i],
          symlink = is_link,
          file = !is_link && file.exists(target) && !dir.exists(target),
          directory = !is_link && dir.exists(target),
          directory_tree = !is_link && dir.exists(target),
          FALSE
        )
        if (!valid_type) {
          stop("Cleanup target is missing or has an unexpected type: ", target,
               call. = FALSE)
        }
      }
      active_index <- NA_integer_

      depth <- lengths(strsplit(artifacts$relative_path, "/", fixed = TRUE))
      type_order <- match(
        artifacts$artifact_type,
        c("file", "symlink", "directory", "directory_tree")
      )
      delete_indices <- delete_indices[order(type_order[delete_indices],
                                             -depth[delete_indices])]
      for (i in delete_indices) {
        active_index <- i
        recursive <- artifacts$artifact_type[i] == "directory_tree"
        status <- unlink(targets[i], recursive = recursive, force = FALSE)
        if (status != 0 || path_exists_without_following(targets[i])) {
          stop("Failed to delete current-run artifact: ", targets[i], call. = FALSE)
        }
        audit$action[i] <- "deleted"
      }
      active_index <- NA_integer_
    }

    write_tsv_atomic(audit, audit_path)
    invisible(audit)
  }, error = function(error) {
    failed_rows <- audit$action == "pending"
    if (!is.na(active_index)) failed_rows[active_index] <- TRUE
    if (!any(failed_rows)) failed_rows <- audit$action != "deleted"
    audit$action[failed_rows] <- "failed"
    audit$reason[failed_rows] <- paste0(
      audit$reason[failed_rows], "; ", conditionMessage(error)
    )
    write_tsv_atomic(audit, audit_path)
    stop(error)
  })
}

# Leaf labels are complete NCBI taxon names with spaces encoded as underscores.
# Legacy labels ending in an assembly accession, including a missing version, are rejected.
ACCESSION_RE <- "(^|_)GC[AF]_[0-9]+(\\.[0-9]+)?$"

# "Chaunax_sp._Z400" -> "Chaunax sp. Z400". Preserve all taxon-name components so
# strains, subspecies, and informal species identifiers remain distinct Stage-0 taxa.
taxon_from_label <- function(name) {
  if (!is.na(name) && grepl("^'.*'$", name)) {
    name <- substr(name, 2, nchar(name) - 1)
    name <- gsub("''", "'", name, fixed = TRUE)
  }
  if (is.na(name) || grepl("[[:space:]]", name) || grepl(ACCESSION_RE, name)) {
    return(NA_character_)
  }
  tokens <- strsplit(name, "_", fixed = TRUE)[[1]]
  if (length(tokens) >= 2 && all(nzchar(tokens))) paste(tokens, collapse = " ")
  else NA_character_
}

species_from_label <- function(labels) {
  vapply(labels, taxon_from_label, character(1), USE.NAMES = FALSE)
}

# Short names must be stable per species and must NOT carry a trio slot number
# (sbst_fromDwl.sh's make_short_name appends the slot). Used only for readable
# output/log naming (selected_trios.tsv's out_short/in1_short/in2_short columns
# etc.) -- the last-train cache is keyed on NCBI accession, not on this name.
make_short_names <- function(species) {
  base <- vapply(species, function(name) {
    tokens <- strsplit(name, "[ _]+")[[1]]
    first <- substr(tokens[1], 1, 3)
    second <- if (length(tokens) >= 2) substr(tokens[2], 1, 3) else substr(tokens[1], 4, 6)
    paste0(first, second)
  }, character(1), USE.NAMES = FALSE)

  # Disambiguate collisions so two species never show the same short name in output.
  for (dup in unique(base[duplicated(base)])) {
    idx <- which(base == dup)
    base[idx] <- paste0(dup, seq_along(idx))
  }
  base
}

descendant_tips <- function(tree, node) {
  # A tip is its own only descendant tip; otherwise recurse into the children.
  # (Without the tip base case a single-leaf sister clade is silently dropped, so
  # a pair whose only outgroups are lone leaves would get no candidates at all.)
  if (node <= Ntip(tree)) return(node)
  kids <- tree$edge[tree$edge[, 1] == node, 2]
  unlist(lapply(kids, descendant_tips, tree = tree))
}

sister_tips <- function(tree, node) {
  parent <- tree$edge[tree$edge[, 2] == node, 1]
  if (!length(parent)) return(integer(0))
  siblings <- setdiff(tree$edge[tree$edge[, 1] == parent, 2], node)
  unlist(lapply(siblings, descendant_tips, tree = tree))
}

# --- Stage 0: prune the tree to one best assembly per species ----------------

# Fetch every current assembly of each species on the tree (not just the leaf accession),
# hard-filter unsuitable candidates, build a metadata baseline, attach QC for audit and an
# optional explicit strict override, and prune the tree to species that keep an assembly.
# The tree is topology only; leaves$accession becomes the chosen accession, which may differ
# from the leaf label's accession (intentional -- a better off-tree assembly, or recovery of
# a species whose leaf accession is dead).  Duplicate-species leaves collapse to one tip.
prune_to_best_assembly <- function(tree, out_dir, min_rel_contig_n50 = 0,
                                   stage0_top_k = 3L, assembly_qc = NULL,
                                   allow_qc_override = FALSE) {
  tip_species <- species_from_label(tree$tip.label)
  invalid_labels <- is.na(tip_species) | !nzchar(tip_species)
  if (any(invalid_labels)) {
    stop(
      "Stage 0 requires accession-free complete taxon names with spaces encoded as underscores. ",
      "Remove trailing NCBI accessions with src/select/strip_newick_accessions.py. ",
      "Invalid label(s): ",
      paste(head(tree$tip.label[invalid_labels], 5), collapse = ", "),
      call. = FALSE
    )
  }
  taxa <- unique(tip_species)

  tax_file <- file.path(out_dir, "tree_taxa.txt")
  meta_file <- file.path(out_dir, "assembly_metadata.tsv")
  writeLines(taxa, tax_file)

  log_stage("Stage 0: fetching all current assemblies for", length(taxa), "species")
  run("python3", c(shQuote(FETCH_METADATA),
                   "--taxons", shQuote(tax_file),
                   "--out", shQuote(meta_file)))

  meta <- read_tsv(meta_file)
  external_qc <- NULL
  if (!is.null(assembly_qc)) {
    if (!file.exists(assembly_qc)) {
      stop("External assembly QC file does not exist: ", assembly_qc, call. = FALSE)
    }
    external_qc <- read_tsv(assembly_qc)
  }
  ranked <- rank_assembly_candidates(
    meta,
    min_rel_contig_n50 = min_rel_contig_n50,
    shortlist_k = stage0_top_k,
    external_qc = external_qc,
    allow_qc_override = allow_qc_override
  )
  best <- ranked$best
  write_tsv(ranked$audit, file.path(out_dir, "assembly_candidates_audit.tsv"))
  write_tsv(ranked$shortlist, file.path(out_dir, "assembly_shortlist.tsv"))
  review_file <- file.path(out_dir, "assembly_review.tsv")
  write_tsv(ranked$review, review_file)
  review_species <- length(unique(ranked$review$species))
  log_stage("  ", review_species, "species require assembly review ->", review_file)

  # Keep one tip per species (first in tree order); duplicate-species tips are dropped.
  leaves <- tibble::tibble(
    tip = tree$tip.label,
    species = tip_species,
    .tip_index = seq_along(tree$tip.label)
  ) %>%
    inner_join(best, by = "species") %>%
    group_by(species) %>%
    slice_head(n = 1) %>%
    ungroup()

  # Report dropped species, distinguishing a name mismatch (NCBI returned nothing under the
  # label-parsed name) from a species whose candidates all failed Stage-0 quality screening.
  fetched_species <- unique(meta$species)
  dropped_species <- setdiff(unique(tip_species), unique(leaves$species))
  no_metadata <- setdiff(dropped_species, fetched_species)
  quality_dropped <- setdiff(dropped_species, no_metadata)
  if (length(no_metadata)) {
    log_stage("  dropping", length(no_metadata),
              "species with NO matching NCBI metadata (name mismatch or no assemblies):",
              paste(no_metadata, collapse = ", "))
  }
  if (length(quality_dropped)) {
    log_stage("  dropping", length(quality_dropped),
              "species with no assembly passing Stage-0 quality screening:",
              paste(quality_dropped, collapse = ", "))
  }

  # Note collapsed duplicate-species leaves: one tip (first in tree order) is kept, so the
  # caller knows a topological position among duplicates was chosen arbitrarily.
  dup_kept <- intersect(unique(tip_species[duplicated(tip_species)]), leaves$species)
  if (length(dup_kept)) {
    log_stage("  note:", length(dup_kept),
              "species on multiple leaves; kept one tip each (first in tree order):",
              paste(dup_kept, collapse = ", "))
  }

  if (!nrow(leaves)) {
    stop("Stage 0 kept 0 species: every leaf lacked matching NCBI metadata or every candidate ",
         "failed quality screening. Check assembly_candidates_audit.tsv and the tree's species ",
         "names.", call. = FALSE)
  }

  gate_note <- if (min_rel_contig_n50 > 0)
    sprintf("(relative contig-N50 gate >= %.4g)", min_rel_contig_n50) else ""
  log_stage("  ", nrow(meta), "assembly records ->", nrow(ranked$shortlist),
            "shortlisted ->", nrow(leaves), "species kept", gate_note)

  pruned <- drop.tip(tree, setdiff(seq_along(tree$tip.label), leaves$.tip_index))
  leaves <- leaves[match(pruned$tip.label, leaves$tip), ]
  leaves$.tip_index <- NULL
  leaves$short_name <- make_short_names(leaves$species)

  list(tree = pruned, leaves = leaves)
}

# --- Stage 1: per ingroup pair, ordered candidate outgroups (near -> far) -----

# For the ingroup couple (tipA, tipB), list candidate outgroups nearest first:
# the sister clade of their MRCA, then the sister of the grandparent, ... up to
# the root.  On an ultrametric tree the members of one sister clade are tied, so
# they are ordered by Stage-0 assembly quality (contig N50) as a tie-break.
# Returns an empty vector when the pair's MRCA is the root (no outgroup exists).
candidate_outgroups <- function(tree, tipA, tipB, leaves) {
  mrca <- getMRCA(tree, c(tipA, tipB))
  root <- Ntip(tree) + 1L
  if (is.null(mrca)) return(character(0))

  cands <- integer(0)
  node <- mrca
  while (length(node) == 1 && node != root) {
    sibs <- sister_tips(tree, node)
    if (length(sibs)) {
      # Relative (size-normalized) contig N50 tie-break: compare contiguity fairly across
      # species with very different genome sizes, not raw base counts.
      q <- leaves$rel_contig_n50[match(tree$tip.label[sibs], leaves$tip)]
      q[is.na(q)] <- 0
      cands <- c(cands, sibs[order(-q)])
    }
    parent <- tree$edge[tree$edge[, 2] == node, 1]
    node <- if (length(parent)) parent else root
  }
  tree$tip.label[cands]
}

# Choose which ingroup couples to consider, using tree topology only (no last-train).
#
#   "matching" (default): a phylogenetically-independent sister-pair design -- pair each
#     species with its closest relative and use each species at most once.  Candidate
#     couples are ranked closest-first by TOPOLOGICAL distance (the number of edges between
#     the tips: every edge counts as 1, so a cherry = distance 2 and beats any deeper pair)
#     and greedily matched, giving ~n/2 disjoint, independent sister-pair couples.  Topology
#     (not branch length) is used deliberately: rate variation or additive (substitutions/
#     site) branch lengths could otherwise rank a non-sister pair with short branches ahead
#     of a true cherry with long branches.  Only couples with an available outgroup (MRCA is
#     not the root) are eligible, so a species is never spent on a couple that could yield no
#     trio.  A species left unmatched on a ladder-like tree simply produces no trio.
#
#   "all": every tip pair (exhaustive C(n,2)); the pre-restriction behavior, kept as an
#     opt-out / validation path.
#
# Returns a list of c(tipA, tipB) character pairs.
ingroup_pairs <- function(tree, leaves, mode = "matching") {
  tips <- tree$tip.label
  n <- length(tips)
  if (n < 2) return(list())

  if (mode == "all") {
    combs <- utils::combn(n, 2)
    return(lapply(seq_len(ncol(combs)), function(k) tips[combs[, k]]))
  }

  # Closeness is measured on the tree TOPOLOGY, not branch lengths: set every edge to 1
  # so a pair's distance is the number of edges between them.  A cherry (true sisters) is
  # always distance 2 -- the minimum -- so cherries are matched before any non-sister pair.
  # Real branch lengths would break this: on a tree with rate variation or additive
  # (substitutions/site) lengths a non-sister pair with short branches can have a smaller
  # patristic distance than a true sister pair with long branches, so branch-length ranking
  # could pair non-sisters.  Topology also makes the matching independent of the tree's units
  # (time vs substitutions) and of whether it is ultrametric.
  t2 <- tree
  t2$edge.length <- rep(1, nrow(t2$edge))
  D <- cophenetic(t2)[tips, tips]
  M <- mrca(tree)
  root <- Ntip(tree) + 1L

  ij <- which(upper.tri(D), arr.ind = TRUE)     # candidate i<j pairs (by tip index)
  eligible <- M[ij] != root                     # keep only couples that have an outgroup
  ij <- ij[eligible, , drop = FALSE]
  if (!nrow(ij)) return(list())

  # Deterministic ordering: fewest edges (closest / most-nested) first, then better combined
  # assembly quality (relative contig N50, size-normalized), then tip names -- so the same
  # tree always yields the same matching.
  q <- leaves$rel_contig_n50[match(tips, leaves$tip)]
  q[is.na(q)] <- 0
  score <- q[ij[, 1]] + q[ij[, 2]]
  ord <- order(D[ij], -score, tips[ij[, 1]], tips[ij[, 2]])

  used <- logical(n)
  pairs <- list()
  for (k in ord) {
    a <- ij[k, 1]; b <- ij[k, 2]
    if (used[a] || used[b]) next
    used[a] <- TRUE; used[b] <- TRUE
    pairs[[length(pairs) + 1L]] <- c(tips[a], tips[b])
  }
  pairs
}

# Enumerate the selected ingroup couples with their ordered candidate-outgroup list.
# Used by --dry-run to expose the search order without training.
enumerate_pairs_dry <- function(tree, leaves, mode = "matching") {
  rows <- list()
  for (p in ingroup_pairs(tree, leaves, mode)) {
    cands <- candidate_outgroups(tree, p[1], p[2], leaves)
    if (!length(cands)) next
    rows[[length(rows) + 1]] <- data.frame(
      in1_tip = p[1], in2_tip = p[2],
      n_candidates = length(cands),
      candidate_order = paste(cands, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Mirror of trio_filter.py pattern5: the outgroup is congeneric with exactly one
# ingroup species.  The thesis rule excludes such a trio whatever the identities
# say, so skipping it here avoids a wasted last-train.
is_two_vs_one <- function(g_out, g_in1, g_in2) {
  o <- tolower(g_out); a <- tolower(g_in1); b <- tolower(g_in2)
  xor(o == a, o == b)
}

# --- on-demand download + last-train, memoized ------------------------------

train_file_path <- function(cache_dir, acc_a, acc_b, date) {
  file.path(cache_dir, sprintf("%s2%s_%s.train", acc_a, acc_b, date))
}

parse_train_identity <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  identity_lines <- which(startsWith(lines, IDENTITY_PREFIX))
  if (!length(identity_lines)) return(NA_real_)

  # Intermediate training iterations also contain identity lines. The final
  # identity is the last one and is valid only when followed by the complete
  # parameter block consumed by lastal -p.
  identity_idx <- tail(identity_lines, 1)
  value_text <- trimws(sub(IDENTITY_PREFIX, "", lines[identity_idx], fixed = TRUE))
  if (!grepl("^[0-9]+([.][0-9]+)?$", value_text)) return(NA_real_)
  value <- as.numeric(value_text)
  if (!is.finite(value) || value < 0 || value > 100) return(NA_real_)
  if (identity_idx >= length(lines)) return(NA_real_)

  final_block <- lines[(identity_idx + 1):length(lines)]
  required_parameters <- c(
    "#last -t", "#last -a", "#last -A",
    "#last -b", "#last -B", "#last -S"
  )
  if (!all(vapply(
    required_parameters,
    function(prefix) any(startsWith(final_block, prefix)),
    logical(1)
  ))) return(NA_real_)

  matrix_headers <- which(
    final_block == "# score matrix (query letters = columns, reference letters = rows):"
  )
  if (!length(matrix_headers)) return(NA_real_)
  matrix_idx <- tail(matrix_headers, 1)
  if (matrix_idx >= length(final_block)) return(NA_real_)

  score_number <- "-?[0-9]+([.][0-9]+)?"
  row_pattern <- paste0(
    "^([ACGT])[[:space:]]+", score_number,
    "[[:space:]]+", score_number,
    "[[:space:]]+", score_number,
    "[[:space:]]+", score_number,
    "[[:space:]]*$"
  )
  matrix_rows <- sub(
    row_pattern, "\\1",
    grep(row_pattern, final_block[(matrix_idx + 1):length(final_block)], value = TRUE)
  )
  if (!setequal(matrix_rows, c("A", "C", "G", "T"))) return(NA_real_)

  value
}

# Build closures that download a genome and last-train a species pair at most
# once each (memoized), so a species trained for one ingroup pair is reused for
# every other.  The shared cache dir also lets last_train.sh skip work across runs.
# Cached/memoized by accession pair, not species, so a stale cache entry can never
# be silently reused after Stage 0 picks a different assembly for the same species
# in a later run.
make_fetchers <- function(leaves, opts) {
  fasta_cache <- new.env(parent = emptyenv())
  idt_cache <- new.env(parent = emptyenv())
  owned_artifacts <- new.env(parent = emptyenv())
  owned_artifacts$records <- empty_owned_artifacts()
  counters <- new.env(parent = emptyenv())
  counters$downloads <- 0L
  counters$trains <- 0L

  cache_dir <- if (!is.null(opts$train_cache_dir) && nzchar(opts$train_cache_dir)) {
    opts$train_cache_dir
  } else {
    file.path(opts$out_dir, "train_cache")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(THREAD_NUM = opts$threads)

  append_owned <- function(records) {
    if (nrow(records)) {
      owned_artifacts$records <- rbind(owned_artifacts$records, records)
    }
  }

  get_fasta <- function(acc) {
    hit <- fasta_cache[[acc]]
    if (!is.null(hit)) return(hit)
    download_args <- c(shQuote(DWL_ORGANISM), acc, "--out-dir", shQuote(opts$genome_dir))
    manifest_path <- NULL
    if (!opts$keep_unused_species_data) {
      manifest_path <- tempfile(pattern = ".trio-artifacts-", fileext = ".tsv")
      download_args <- c(
        download_args, "--artifact-manifest", shQuote(manifest_path)
      )
      on.exit(unlink(manifest_path), add = TRUE)
    }
    result <- run("bash", download_args)
    if (!is.null(manifest_path)) {
      append_owned(read_artifact_manifest(manifest_path, acc))
    }
    # dwl_organism.sh prints dir_name|fasta_path|summary_json|raw_name|ncbi_full|tax_json
    fields <- strsplit(tail(result, 1), "|", fixed = TRUE)[[1]]
    if (length(fields) < 2 || !nzchar(fields[2])) {
      stop("dwl_organism.sh returned no FASTA path for ", acc, call. = FALSE)
    }
    # last_train.sh cds into its lastdb directory, so the FASTA path must be absolute.
    path <- normalizePath(fields[2], mustWork = TRUE)
    fasta_cache[[acc]] <- path
    counters$downloads <- counters$downloads + 1L
    path
  }

  get_identity <- function(acc_a, acc_b) {
    # Canonical (unordered) orientation, keyed on NCBI accession -- not species
    # short name -- so a stale cache entry can never be reused for a species
    # whose Stage 0 pick changed to a different assembly between runs. An
    # accession pair uniquely identifies exactly the two genome sequences that
    # were last-trained.
    if (acc_a <= acc_b) {
      a1 <- acc_a; a2 <- acc_b
    } else {
      a1 <- acc_b; a2 <- acc_a
    }
    key <- paste(a1, a2, sep = "|")
    hit <- idt_cache[[key]]
    if (!is.null(hit)) return(hit)

    train <- train_file_path(cache_dir, a1, a2, opts$date)
    val <- parse_train_identity(train)
    if (is.na(val)) {
      # last_train.sh treats every existing path as complete. Remove partial
      # files and dangling links first so interrupted cache writes can recover.
      link_target <- Sys.readlink(train)
      path_exists <- file.exists(train) || (!is.na(link_target) && nzchar(link_target))
      if (path_exists) {
        unlink_status <- unlink(train)
        link_target <- Sys.readlink(train)
        path_exists <- file.exists(train) || (!is.na(link_target) && nzchar(link_target))
        if (unlink_status != 0 || path_exists) {
          stop("Could not remove invalid last-train cache file: ", train, call. = FALSE)
        }
      }
      fa1 <- get_fasta(a1); fa2 <- get_fasta(a2)
      db_relative <- paste0(a1, "db_", opts$date)
      db_path <- file.path(cache_dir, db_relative)
      db_preexisting <- path_exists_without_following(db_path)
      run("bash", c(shQuote(LAST_TRAIN), opts$date, shQuote(fa1), shQuote(fa2),
                    a1, a2, shQuote(cache_dir)))
      counters$trains <- counters$trains + 1L
      val <- parse_train_identity(train)
      if (!opts$keep_unused_species_data) {
        train_link <- Sys.readlink(train)
        train_is_link <- !is.na(train_link) && nzchar(train_link)
        if (file.exists(train) && !dir.exists(train) && !train_is_link) {
          append_owned(data.frame(
            accession_1 = a1,
            accession_2 = a2,
            storage = "train_cache",
            artifact_type = "file",
            relative_path = basename(train),
            stringsAsFactors = FALSE
          ))
        }
        db_link <- Sys.readlink(db_path)
        db_is_link <- !is.na(db_link) && nzchar(db_link)
        if (!db_preexisting && dir.exists(db_path) && !db_is_link) {
          append_owned(data.frame(
            accession_1 = a1,
            accession_2 = "",
            storage = "train_cache",
            artifact_type = "directory_tree",
            relative_path = db_relative,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    idt_cache[[key]] <- val
    val
  }

  list(
    get_identity = get_identity,
    counters = counters,
    get_owned_artifacts = function() owned_artifacts$records
  )
}

# --- the thesis rule, judged one candidate at a time ------------------------

# Shell one candidate trio through trio_filter.py (the single source of the
# thesis rule) and return its verdict row.
evaluate_one <- function(row, opts, tmp_in, tmp_out) {
  write_tsv(row, tmp_in)
  run("python3", c(shQuote(TRIO_FILTER),
                   "--trios", shQuote(tmp_in),
                   "--idt-threshold", opts$idt_threshold,
                   "--out", shQuote(tmp_out)),
      quiet_stderr = TRUE)
  read_tsv(tmp_out)[1, , drop = FALSE]
}

build_row <- function(infoO, infoA, infoB, idt_12, idt_13, idt_23, tries) {
  data.frame(
    out_tip = infoO$tip, in1_tip = infoA$tip, in2_tip = infoB$tip,
    out_acc = infoO$accession, in1_acc = infoA$accession, in2_acc = infoB$accession,
    out_species = infoO$species, in1_species = infoA$species, in2_species = infoB$species,
    out_short = infoO$short_name, in1_short = infoA$short_name, in2_short = infoB$short_name,
    # Slot 1 is the outgroup, matching org1/org2/org3 in the rest of the pipeline.
    genus_1 = infoO$genus, genus_2 = infoA$genus, genus_3 = infoB$genus,
    idt_12 = idt_12, idt_13 = idt_13, idt_23 = idt_23,
    outgroup_tries = tries,
    stringsAsFactors = FALSE
  )
}

# --- Stage 1 driver: one nearest-but-external outgroup per ingroup pair -------

select_trios <- function(tree, leaves, opts, fetch = NULL) {
  if (is.null(fetch)) fetch <- make_fetchers(leaves, opts)
  get_identity <- fetch$get_identity
  threshold <- opts$idt_threshold
  tmp_in <- file.path(opts$out_dir, ".candidate_in.tsv")
  tmp_out <- file.path(opts$out_dir, ".candidate_out.tsv")

  info_of <- function(tip) leaves[match(tip, leaves$tip), , drop = FALSE]

  couples <- ingroup_pairs(tree, leaves, opts$ingroup_pairing)
  covered <- unique(unlist(couples))
  uncovered <- setdiff(tree$tip.label, covered)
  results <- list()
  no_outgroup <- 0L
  too_diverged <- 0L
  no_valid_outgroup <- 0L

  log_stage(sprintf("Stage 1: %s pairing -> %d ingroup couple(s); %d species not paired",
                    opts$ingroup_pairing, length(couples), length(uncovered)))
  log_stage("Stage 1: searching one outgroup per ingroup pair (nearest-first, first-pass)")

  for (couple in couples) {
    tipA <- couple[1]; tipB <- couple[2]
    cands <- candidate_outgroups(tree, tipA, tipB, leaves)
    if (!length(cands)) { no_outgroup <- no_outgroup + 1L; next }

    infoA <- info_of(tipA); infoB <- info_of(tipB)
    idt_23 <- get_identity(infoA$accession, infoB$accession)
    if (is.na(idt_23) || idt_23 <= threshold) { too_diverged <- too_diverged + 1L; next }

    chosen <- NULL
    tries <- 0L
    for (og in cands) {
      infoO <- info_of(og)
      # Free pre-filter: a two-vs-one genus split is always excluded -> never train it.
      if (is_two_vs_one(infoO$genus, infoA$genus, infoB$genus)) next
      if (tries >= opts$max_outgroup_tries) break
      tries <- tries + 1L

      idt_12 <- get_identity(infoO$accession, infoA$accession)
      idt_13 <- get_identity(infoO$accession, infoB$accession)
      if (is.na(idt_12) || is.na(idt_13)) next  # unusable alignment, try the next candidate

      row <- build_row(infoO, infoA, infoB, idt_12, idt_13, idt_23, tries)
      verdict <- evaluate_one(row, opts, tmp_in, tmp_out)

      if (identical(verdict$excluded, "FALSE")) {           # pass -> fix this trio, stop
        chosen <- verdict
        break
      }
      if (identical(verdict$idt_threshold_condition, "FALSE")) {
        # too far: an outgroup-ingroup identity fell below the threshold, and every
        # farther candidate is more diverged still -> the window has closed.
        break
      }
      # else "too close" (strict ordering failed) -> the next candidate is farther out.
    }

    if (!is.null(chosen)) {
      results[[length(results) + 1]] <- chosen
      log_stage(sprintf("  %s + %s  ->  %s  (idt 12/13/23 = %s/%s/%s; %d tried)",
                        infoA$species, infoB$species, chosen$out_species,
                        chosen$idt_12, chosen$idt_13, chosen$idt_23, tries))
    } else {
      no_valid_outgroup <- no_valid_outgroup + 1L
    }
  }

  log_stage(sprintf("last-train runs: %d | downloads: %d", fetch$counters$trains,
                    fetch$counters$downloads))
  log_stage(sprintf("ingroup couples: %d selected | %d no outgroup on tree | %d couple below threshold | %d no external outgroup found | %d species unpaired",
                    length(results), no_outgroup, too_diverged, no_valid_outgroup,
                    length(uncovered)))

  if (!length(results)) return(data.frame())
  do.call(rbind, results)
}

# --- main --------------------------------------------------------------------

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  load_runtime_dependencies()
  dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

  tree <- read.tree(opts$tree)
  log_stage("tree:", opts$tree, "-", Ntip(tree), "leaves")

  pruned <- prune_to_best_assembly(
    tree,
    opts$out_dir,
    opts$min_rel_contig_n50,
    opts$stage0_top_k,
    opts$assembly_qc,
    opts$allow_qc_override
  )
  write_tsv(pruned$leaves, file.path(opts$out_dir, "selected_assemblies.tsv"))
  tree <- pruned$tree
  leaves <- pruned$leaves

  if (opts$dry_run) {
    pairs <- enumerate_pairs_dry(tree, leaves, opts$ingroup_pairing)
    out <- file.path(opts$out_dir, "candidate_outgroups.tsv")
    write_tsv(pairs, out)
    log_stage("Stage 1 (dry-run):", nrow(pairs), "ingroup couple(s) with an available",
              "outgroup (", opts$ingroup_pairing, "pairing); candidate order ->", out)
    return(invisible(NULL))
  }

  fetch <- make_fetchers(leaves, opts)
  selected <- select_trios(tree, leaves, opts, fetch)
  selected_file <- file.path(opts$out_dir, "selected_trios.tsv")
  write_tsv(selected, selected_file)
  log_stage("Selected", nrow(selected), "trios (one per ingroup pair) ->", selected_file)

  if (!opts$keep_unused_species_data) {
    cache_dir <- if (!is.null(opts$train_cache_dir) && nzchar(opts$train_cache_dir)) {
      opts$train_cache_dir
    } else {
      file.path(opts$out_dir, "train_cache")
    }
    cleanup_file <- file.path(opts$out_dir, "cleanup_unused_species.tsv")
    cleanup_owned_artifacts(
      fetch$get_owned_artifacts(), selected, opts$genome_dir, cache_dir, cleanup_file
    )
    log_stage("Current-run unused-species cleanup audit ->", cleanup_file)
  }
}

if (sys.nframe() == 0) {
  main()
}
