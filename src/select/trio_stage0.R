# Helpers loaded by trio_selection.R into its calling environment.

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
