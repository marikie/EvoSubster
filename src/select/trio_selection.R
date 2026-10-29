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
# thesis rule. Both ingroup-minus-outgroup identity gaps must exceed
# --min-idt-gap-between-out-in, unless the target genus pattern exempts them.
# Separate ingroup/outgroup identity thresholds set the "not too far" boundary.
# Nearest-first + first-pass trains as few pairs as possible.
#
# Stage 0  prune the tree to one assembly per species      [NCBI metadata, no download]
#          write pruned_tree.nwk/.pdf with accession labels (also for --dry-run)
# Stage 1  per ingroup pair: search outgroups near->far,   [dwl_organism.sh + last_train.sh
#          stop at the first thesis-rule pass -> one trio    on demand; trio_filter.py to judge]
#          write selected_trios.nwk/.pdf with trio roles when trios are selected (stale ones from a previous run are removed)
#
# Usage:
#   Rscript src/select/trio_selection.R --tree Primates_tree.nwk [options]
#
# Depends on: ape, dplyr (R); R built with cairo (capabilities("cairo")); python3, LAST, curl, jq (external).

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
source(ASSEMBLY_SELECTION, local = environment())
source(file.path(SCRIPT_DIR, "selection_execution.R"), local = environment())
source(file.path(SCRIPT_DIR, "trio_options.R"), local = environment())
source(file.path(SCRIPT_DIR, "trio_io.R"), local = environment())
source(file.path(SCRIPT_DIR, "trio_cleanup.R"), local = environment())
source(file.path(SCRIPT_DIR, "trio_stage0.R"), local = environment())
source(file.path(SCRIPT_DIR, "trio_search.R"), local = environment())
source(file.path(SCRIPT_DIR, "trio_cache.R"), local = environment())
source(file.path(SCRIPT_DIR, "trio_tree_outputs.R"), local = environment())

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

# --- main --------------------------------------------------------------------

main <- function(argv = commandArgs(trailingOnly = TRUE), script_path = NULL) {
  opts <- parse_args(argv)
  if (is.null(script_path)) {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) script_path <- sub("^--file=", "", file_arg[1])
  }
  execution_info <- record_selection_execution(opts, argv, script_path)
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
  stage0_files <- write_stage0_tree(tree, leaves, opts$out_dir)
  log_stage("Stage 0 tree ->", stage0_files[["nwk"]], "/", stage0_files[["pdf"]])
  # Stage 0 output is current from here on, so Stage 1 tree files from a previous run are
  # stale whether this run stops at --dry-run or selects zero trios: clear them now and
  # write them again only when trios are selected.
  remove_stage1_tree_outputs(opts$out_dir)

  if (opts$dry_run) {
    pairs <- enumerate_pairs_dry(tree, leaves, opts$ingroup_pairing)
    out <- file.path(opts$out_dir, "candidate_outgroups.tsv")
    write_tsv(pairs, out)
    log_stage("Stage 1 (dry-run):", nrow(pairs), "ingroup couple(s) with an available",
              "outgroup (", opts$ingroup_pairing, "pairing); candidate order ->",
              paste0(out, "; no Stage 1 tree files"))
    return(invisible(NULL))
  }

  fetch <- make_fetchers(leaves, opts)
  selected <- select_trios(tree, leaves, opts, fetch)
  selected_file <- file.path(opts$out_dir, "selected_trios.tsv")
  write_tsv(selected, selected_file)
  if (nrow(selected)) {
    stage1_files <- write_selected_trio_tree(tree, leaves, selected, opts$out_dir, execution_info)
    log_stage("Stage 1 tree ->", stage1_files[["nwk"]], "/", stage1_files[["pdf"]])
  } else {
    log_stage("Stage 1 tree: no trios selected; none written")
  }
  log_stage("Selected", nrow(selected), "trios (one per ingroup pair) ->", selected_file)

  if (!opts$keep_unused_species_data) {
    cache_dir <- resolve_train_cache_dir(opts)
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
