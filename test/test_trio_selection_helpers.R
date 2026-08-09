#!/usr/bin/env Rscript

# Unit tests for dependency-free trio-selection helpers. The package loader is
# temporarily suppressed so these parsing/path tests do not require ape or dplyr.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))

real_library <- library
library <- function(...) invisible(TRUE)
source(file.path(root, "src", "select", "trio_selection.R"))
library <- real_library

fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) {
    cat("ok:", msg, "\n")
  } else {
    cat("FAIL:", msg, "\n")
    fail <<- fail + 1L
  }
}

labels <- c(
  "Homo_sapiens",
  "Chaunax_sp._Z400",
  "'Chaunax_sp._Z400'",
  "Pseudoliparis_sp._2_HX-2024",
  "Gasterosteus_aculeatus_aculeatus",
  "Drosophila_melanogaster_strain_A"
)
check(
  identical(
    unname(species_from_label(labels)),
    c(
      "Homo sapiens",
      "Chaunax sp. Z400",
      "Chaunax sp. Z400",
      "Pseudoliparis sp. 2 HX-2024",
      "Gasterosteus aculeatus aculeatus",
      "Drosophila melanogaster strain A"
    )
  ),
  "species parser preserves every component of accession-free taxon labels"
)

check(
  identical(taxon_from_label("'Chaunax_sp._Z400'"), "Chaunax sp. Z400") &&
    is.na(taxon_from_label("Chaunax_sp._Z400_GCA_037577475")),
  "scalar taxon parser preserves accepted and rejected label behavior"
)

malformed <- species_from_label(c(
  "Homo",
  "_",
  "Homo sapiens",
  "Mus_musculus_GCF_000001635.27",
  "Mus_musculus_GCA_000001635",
  "GCA_000001405.40"
))
check(
  all(is.na(malformed)),
  "species parser rejects malformed, space-containing, and accession-suffixed labels"
)

check(
  identical(
    train_file_path(
      "/tmp/cache",
      "GCA_000000001.1",
      "GCF_000000002.2",
      "20260723"
    ),
    "/tmp/cache/GCA_000000001.12GCF_000000002.2_20260723.train"
  ),
  "train cache filename is keyed by the complete accession pair"
)

help_text <- paste(capture.output(usage()), collapse = "\n")
check(
  grepl("accession-free complete taxon names", help_text, fixed = TRUE) &&
    grepl("strip_newick_accessions.py", help_text, fixed = TRUE),
  "tree help documents the canonical complete taxon labels and legacy converter"
)
check(
  grepl("--out-dir DIR[[:space:]]+Output directory", help_text),
  "--out-dir help does not imply that every cache must live there"
)
check(
  grepl("default: <out-dir>/train_cache", help_text, fixed = TRUE),
  "--train-cache-dir help states its default"
)
check(
  grepl("--keep-unused-species-data", help_text, fixed = TRUE),
  "CLI help documents the unused-species cleanup opt-out"
)

help_process <- suppressWarnings(system2(
  file.path(R.home("bin"), "Rscript"),
  c(file.path(root, "src", "select", "trio_selection.R"), "--help"),
  stdout = TRUE,
  stderr = TRUE
))
check(
  is.null(attr(help_process, "status")) && any(grepl("^Usage:", help_process)),
  "CLI help exits successfully without loading runtime packages"
)

stage0_args <- try(
  parse_args(c(
    "--tree", "example.nwk",
    "--stage0-top-k", "4",
    "--assembly-qc", "external_qc.tsv"
  )),
  silent = TRUE
)
check(
  !inherits(stage0_args, "try-error") && identical(stage0_args$stage0_top_k, 4L),
  "Stage 0 CLI accepts a positive per-species shortlist size"
)
check(
  !inherits(stage0_args, "try-error") &&
    identical(stage0_args$assembly_qc, "external_qc.tsv"),
  "Stage 0 CLI accepts an optional external genome-QC table"
)
check(
  grepl("BUSCO genome-mode", help_text, fixed = TRUE),
  "Stage 0 help distinguishes external BUSCO genome-mode from NCBI annotation BUSCO"
)
check(
  grepl("--allow-qc-override", help_text, fixed = TRUE),
  "Stage 0 help documents the explicit strict QC override flag"
)
check(
  grepl(
    "metadata shortlist for QC audit and optional explicit strict override",
    help_text,
    fixed = TRUE
  ),
  "Stage 0 help describes the audit/optional-override shortlist policy"
)
check(
  identical(defaults$allow_qc_override, FALSE),
  "Stage 0 QC override is disabled by default"
)

legacy_tree_errors <- lapply(
  c(
    "Chaunax_sp._Z400_GCA_037577475.1",
    "Chaunax_sp._Z400_GCA_037577475"
  ),
  function(label) try(
    prune_to_best_assembly(
      list(tip.label = label),
      tempfile("trio-stage0-invalid-label-")
    ),
    silent = TRUE
  )
)
check(
  all(vapply(legacy_tree_errors, function(error) {
    inherits(error, "try-error") &&
      grepl("strip_newick_accessions.py", as.character(error), fixed = TRUE)
  }, logical(1))),
  "Stage 0 rejects versioned and versionless accession suffixes with converter guidance"
)

override_args <- try(
  parse_args(c("--tree", "example.nwk", "--allow-qc-override")),
  silent = TRUE
)
check(
  !inherits(override_args, "try-error") &&
    identical(override_args$allow_qc_override, TRUE),
  "Stage 0 CLI accepts the explicit strict QC override flag"
)

# Exercise the prune-to-ranker integration boundary without loading phylogeny
# packages or contacting NCBI. The ranker sentinel distinguishes a forwarded
# TRUE value from an omitted/default FALSE value.
stage0_out <- tempfile("trio-stage0-forwarding-")
dir.create(stage0_out, recursive = TRUE)
write.table(
  data.frame(species = "Homo sapiens", accession = "GCF_000001405.40"),
  file.path(stage0_out, "assembly_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
fake_tree <- list(tip.label = "Homo_sapiens")
real_run <- run
real_rank_assembly_candidates <- rank_assembly_candidates
run <- function(...) character()
rank_assembly_candidates <- function(..., allow_qc_override = FALSE) {
  stop(
    paste0("ranker_allow_qc_override=", allow_qc_override),
    call. = FALSE
  )
}
forwarding_err <- try(
  prune_to_best_assembly(
    fake_tree,
    stage0_out,
    allow_qc_override = TRUE
  ),
  silent = TRUE
)
check(
  inherits(forwarding_err, "try-error") &&
    grepl("ranker_allow_qc_override=TRUE", as.character(forwarding_err), fixed = TRUE),
  "Stage 0 pruning forwards the explicit QC override to assembly ranking"
)
rank_assembly_candidates <- real_rank_assembly_candidates

# Exercise review output with the real ranker. The synthetic metadata mirrors
# the public fetch TSV order so the expected header is independent of ranker internals.
stage0_metadata_row <- function(accession, reference = FALSE, annotated = FALSE,
                                contig_n50 = 1e6) {
  data.frame(
    accession = accession,
    current_accession = accession,
    organism_name = "Homo sapiens",
    organism_tax_id = "9606",
    species = "Homo sapiens",
    genus = "Homo",
    source_database = if (reference) "SOURCE_DATABASE_REFSEQ" else
      "SOURCE_DATABASE_GENBANK",
    refseq_category = if (reference) "reference genome" else "",
    assembly_level = "Chromosome",
    assembly_type = "haploid",
    release_date = "2025-01-01",
    sequencing_tech = "PacBio",
    assembly_method = "ExampleAssembler",
    is_atypical = "false",
    atypical_warnings = "",
    ani_check_status = "",
    checkm_completeness = NA_real_,
    checkm_contamination = NA_real_,
    busco_lineage = "",
    busco_version = "",
    busco_complete = NA_real_,
    busco_duplicated = NA_real_,
    busco_fragmented = NA_real_,
    busco_missing = NA_real_,
    total_sequence_length = 1e8,
    total_ungapped_length = 0.99e8,
    number_of_contigs = 10,
    contig_n50 = contig_n50,
    number_of_scaffolds = 10,
    scaffold_n50 = contig_n50,
    has_annotation = if (annotated) "true" else "false",
    annotation_provider = if (annotated) "NCBI RefSeq" else "",
    annotation_release_date = if (annotated) "2025-02-01" else "",
    assembly_status = "current",
    paired_accession = "",
    stringsAsFactors = FALSE
  )
}

expected_review_columns <- c(
  "accession", "current_accession", "organism_name", "organism_tax_id",
  "species", "genus", "source_database", "refseq_category", "assembly_level",
  "assembly_type", "release_date", "sequencing_tech", "assembly_method",
  "is_atypical", "atypical_warnings", "ani_check_status",
  "checkm_completeness", "checkm_contamination", "busco_lineage",
  "busco_version", "busco_complete", "busco_duplicated", "busco_fragmented",
  "busco_missing", "total_sequence_length", "total_ungapped_length",
  "number_of_contigs", "contig_n50", "number_of_scaffolds", "scaffold_n50",
  "has_annotation", "annotation_provider", "annotation_release_date",
  "assembly_status", "paired_accession", "qc_busco_mode", "qc_busco_lineage",
  "qc_busco_version", "qc_busco_complete", "qc_busco_single",
  "qc_busco_duplicated", "qc_busco_fragmented", "qc_busco_missing",
  "qc_busco_internal_stop", "merqury_qv", "merqury_completeness",
  "rel_contig_n50", "gap_fraction", "is_reference", "is_refseq",
  "assembly_level_rank", "primary_type_rank", "eligible", "exclusion_reason",
  "selection_profile", "shortlist_rank", "shortlisted", "baseline_rank",
  "final_rank", "selected", "selection_basis", "baseline_selected",
  "qc_preferred", "review_required", "review_reason", "override_applied",
  "override_block_reason", "assembly_equivalence_key"
)

write_stage0_fixture <- function(path, metadata) {
  dir.create(path, recursive = TRUE)
  write.table(
    metadata,
    file.path(path, "assembly_metadata.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
})
duplicate_tip_out <- tempfile("trio-stage0-duplicate-tip-")
duplicate_tip_metadata <- rbind(
  transform(
    stage0_metadata_row("GCF_000000101.1", reference = TRUE, annotated = TRUE),
    organism_name = "Genus alpha",
    organism_tax_id = "101",
    species = "Genus alpha",
    genus = "Genus"
  ),
  transform(
    stage0_metadata_row("GCF_000000102.1", reference = TRUE, annotated = TRUE),
    organism_name = "Genus beta",
    organism_tax_id = "102",
    species = "Genus beta",
    genus = "Genus"
  )
)
write_stage0_fixture(duplicate_tip_out, duplicate_tip_metadata)
duplicate_tip_tree <- read.tree(text = "((Genus_alpha:1,Genus_alpha:1):1,Genus_beta:2);")
duplicate_tip_result <- prune_to_best_assembly(duplicate_tip_tree, duplicate_tip_out)
check(
  Ntip(duplicate_tip_result$tree) == 2L &&
    nrow(duplicate_tip_result$leaves) == 2L &&
    !anyDuplicated(duplicate_tip_result$leaves$species),
  "Stage 0 collapses duplicate taxon labels to one tree tip per species"
)

no_review_out <- tempfile("trio-stage0-empty-review-")
no_review_metadata <- rbind(
  stage0_metadata_row("GCF_000001405.40", reference = TRUE, annotated = TRUE),
  stage0_metadata_row("GCA_000001406.1", contig_n50 = 2e6)
)
write_stage0_fixture(no_review_out, no_review_metadata)
invisible(try(prune_to_best_assembly(fake_tree, no_review_out), silent = TRUE))
no_review_file <- file.path(no_review_out, "assembly_review.tsv")
no_review_table <- read.delim(
  no_review_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)
check(
  identical(names(no_review_table), expected_review_columns) &&
    nrow(no_review_table) == 0L,
  "Stage 0 real ranker writes the complete review schema when no review is required"
)

review_out <- tempfile("trio-stage0-populated-review-")
review_metadata <- rbind(
  stage0_metadata_row("GCF_000001405.40", reference = TRUE, annotated = TRUE),
  stage0_metadata_row("GCA_000001406.1", contig_n50 = 2e6)
)
write_stage0_fixture(review_out, review_metadata)
review_qc_file <- file.path(review_out, "external_qc.tsv")
write.table(
  data.frame(
    accession = c("GCF_000001405.40", "GCA_000001406.1"),
    qc_busco_mode = c("genome", "genome"),
    qc_busco_lineage = c("primates_odb10", "primates_odb10"),
    qc_busco_version = c("5.8.2", "5.8.2"),
    qc_busco_complete = c(95, 98),
    qc_busco_single = c(93, 97),
    qc_busco_duplicated = c(2, 1),
    qc_busco_fragmented = c(2, 1),
    qc_busco_missing = c(3, 1),
    stringsAsFactors = FALSE
  ),
  review_qc_file,
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
invisible(try(prune_to_best_assembly(
  fake_tree, review_out, stage0_top_k = 2L, assembly_qc = review_qc_file
), silent = TRUE))
review_file <- file.path(review_out, "assembly_review.tsv")
review_table <- read.delim(
  review_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)
check(
  identical(names(review_table), expected_review_columns) &&
    identical(
      sort(review_table$accession),
      c("GCA_000001406.1", "GCF_000001405.40")
    ),
  "Stage 0 real ranker writes every candidate row for a reviewed species"
)
run <- real_run

# An interrupted cache file must be removed and trained again instead of being
# parsed as NA forever. Stub only the external process boundary; exercise the
# real make_fetchers() cache logic.
cache_root <- tempfile("trio-cache-")
dir.create(cache_root, recursive = TRUE)
fake_fasta <- file.path(cache_root, "genome.fna")
writeLines(c(">seq", "ACGT"), fake_fasta)
acc_a <- "GCA_000000001.1"
acc_b <- "GCA_000000002.1"
cache_file <- train_file_path(cache_root, acc_a, acc_b, "20260723")
writeLines(
  c(
    "# substitution percent identity: 91.0",
    "# count matrix (intermediate training iteration)",
    "# score matrix (query letters = columns, reference letters = rows):"
  ),
  cache_file
)

write_complete_train <- function(path, final_identity) {
  writeLines(
    c(
      "# substitution percent identity: 91.0",
      "# count matrix (intermediate training iteration)",
      paste("# substitution percent identity:", final_identity),
      "# ref letter %: 25.0 25.0 25.0 25.0",
      "# qry letter %: 25.0 25.0 25.0 25.0",
      "#last -t4.5",
      "#last -a 20",
      "#last -A 20",
      "#last -b 1",
      "#last -B 1",
      "#last -S 1",
      "# score matrix (query letters = columns, reference letters = rows):",
      "       A      C      G      T",
      "A      6    -10     -7    -11",
      "C    -10      6    -11     -7",
      "G     -7    -11      6    -10",
      "T    -11     -7    -10      6"
    ),
    path
  )
}

real_run <- run
run <- function(command, args, quiet_stderr = FALSE) {
  if (grepl("dwl_organism[.]sh", args[1])) {
    manifest_idx <- match("--artifact-manifest", args)
    if (!is.na(manifest_idx)) {
      manifest_path <- gsub("^'|'$", "", args[manifest_idx + 1])
      writeLines(
        c(
          "accession\tartifact_type\trelative_path",
          paste(args[2], "file", basename(fake_fasta), sep = "\t")
        ),
        manifest_path
      )
    }
    return(paste("Genus_species", fake_fasta, "", "", "", "", sep = "|"))
  }
  if (grepl("last_train[.]sh", args[1])) {
    dir.create(file.path(cache_root, paste0(acc_a, "db_20260723")))
    write_complete_train(cache_file, 88.5)
    return(character())
  }
  stop("Unexpected command in test stub: ", command, " ", paste(args, collapse = " "))
}

fetch <- make_fetchers(
  data.frame(),
  modifyList(
    defaults,
    list(
      out_dir = tempfile("trio-out-"),
      train_cache_dir = cache_root,
      genome_dir = cache_root,
      date = "20260723"
    )
  )
)
identity <- fetch$get_identity(acc_a, acc_b)
run <- real_run

check(
  identical(identity, 88.5),
  "invalid cached train output is replaced and reparsed"
)
check(
  identical(fetch$counters$trains, 1L),
  "invalid cached train output triggers exactly one fresh training run"
)
owned_after_replacement <- fetch$get_owned_artifacts()
check(
  any(
    owned_after_replacement$storage == "train_cache" &
      owned_after_replacement$artifact_type == "file" &
      owned_after_replacement$relative_path == basename(cache_file)
  ),
  "a regenerated complete train cache is owned by the current run"
)
check(
  any(
    owned_after_replacement$storage == "train_cache" &
      owned_after_replacement$artifact_type == "directory_tree" &
      owned_after_replacement$accession_1 == acc_a &
      owned_after_replacement$accession_2 == "" &
      owned_after_replacement$relative_path == paste0(acc_a, "db_20260723")
  ),
  "a current-run LAST database directory is owned by its reference accession"
)
check(
  any(
    owned_after_replacement$storage == "genome" &
      owned_after_replacement$accession_1 == acc_a &
      owned_after_replacement$relative_path == basename(fake_fasta)
  ),
  "downloader manifest rows are aggregated by accession"
)

complete_file <- tempfile("complete-train-")
write_complete_train(complete_file, 87.25)
check(
  identical(parse_train_identity(complete_file), 87.25),
  "complete train parsing uses the final identity instead of an intermediate iteration"
)

out_of_range_file <- tempfile("out-of-range-train-")
write_complete_train(out_of_range_file, 101)
check(
  is.na(parse_train_identity(out_of_range_file)),
  "train parsing rejects a final identity outside 0-100"
)

# Cleanup targets are synthetic current-run ownership records. Unlisted paths
# represent pre-existing data and must survive regardless of selection.
artifact_rows <- function(...) {
  rows <- list(...)
  do.call(rbind, lapply(rows, function(row) {
    data.frame(
      accession_1 = row[[1]], accession_2 = row[[2]], storage = row[[3]],
      artifact_type = row[[4]], relative_path = row[[5]],
      stringsAsFactors = FALSE
    )
  }))
}

cleanup_root <- tempfile("cleanup-root-")
genome_root <- file.path(cleanup_root, "genomes")
train_root <- file.path(cleanup_root, "train")
dir.create(genome_root, recursive = TRUE)
dir.create(train_root, recursive = TRUE)

dir.create(file.path(genome_root, "selected"))
writeLines("selected", file.path(genome_root, "selected", "genome.fna"))
dir.create(file.path(genome_root, "unused"))
writeLines("unused", file.path(genome_root, "unused", "genome.fna"))
writeLines("pre-existing", file.path(genome_root, "pre-existing.fna"))
writeLines("selected cache", file.path(train_root, "selected.train"))
writeLines("unused cache", file.path(train_root, "unused.train"))
dir.create(file.path(train_root, "unused-db"))
writeLines("index", file.path(train_root, "unused-db", "index.bck"))

owned <- artifact_rows(
  list("GCA_SELECTED_OUT.1", "", "genome", "directory_tree", "selected"),
  list("GCA_UNUSED.1", "", "genome", "directory_tree", "unused"),
  list("GCA_SELECTED_OUT.1", "GCA_SELECTED_IN1.1", "train_cache", "file", "selected.train"),
  list("GCA_SELECTED_OUT.1", "GCA_UNUSED.1", "train_cache", "file", "unused.train"),
  list("GCA_UNUSED.1", "", "train_cache", "directory_tree", "unused-db")
)
selected <- data.frame(
  out_acc = "GCA_SELECTED_OUT.1",
  in1_acc = "GCA_SELECTED_IN1.1",
  in2_acc = "GCA_SELECTED_IN2.1",
  stringsAsFactors = FALSE
)
cleanup_audit <- file.path(cleanup_root, "cleanup.tsv")
cleanup_owned_artifacts(owned, selected, genome_root, train_root, cleanup_audit)

check(
  dir.exists(file.path(genome_root, "selected")),
  "selected accession artifacts are retained"
)
check(
  !file.exists(file.path(genome_root, "unused")),
  "unused current-run genome artifacts are deleted"
)
check(
  file.exists(file.path(genome_root, "pre-existing.fna")),
  "untracked pre-existing genome data are retained"
)
check(
  file.exists(file.path(train_root, "selected.train")) &&
    !file.exists(file.path(train_root, "unused.train")) &&
    !dir.exists(file.path(train_root, "unused-db")),
  "generated train files and LAST databases are removed for unused accessions"
)
cleanup_log <- read.delim(cleanup_audit, sep = "\t", stringsAsFactors = FALSE)
check(
  identical(
    names(cleanup_log),
    c("accessions", "artifact_type", "path", "action", "reason")
  ) && setequal(cleanup_log$action, c("retained", "deleted")),
  "cleanup audit records accessions, artifact type, path, action, and reason"
)

zero_root <- tempfile("cleanup-zero-")
zero_genomes <- file.path(zero_root, "genomes")
zero_train <- file.path(zero_root, "train")
dir.create(file.path(zero_genomes, "candidate"), recursive = TRUE)
dir.create(zero_train, recursive = TRUE)
writeLines("candidate", file.path(zero_genomes, "candidate", "genome.fna"))
writeLines("cache", file.path(zero_train, "candidate.train"))
zero_owned <- artifact_rows(
  list("GCA_ZERO_A.1", "", "genome", "directory_tree", "candidate"),
  list("GCA_ZERO_A.1", "GCA_ZERO_B.1", "train_cache", "file", "candidate.train")
)
cleanup_owned_artifacts(
  zero_owned, data.frame(), zero_genomes, zero_train,
  file.path(zero_root, "cleanup.tsv")
)
check(
  !file.exists(file.path(zero_genomes, "candidate")) &&
    !file.exists(file.path(zero_train, "candidate.train")),
  "zero selected trios delete every current-run-owned artifact"
)

unsafe_root <- tempfile("cleanup-unsafe-")
unsafe_genomes <- file.path(unsafe_root, "genomes")
unsafe_train <- file.path(unsafe_root, "train")
dir.create(unsafe_genomes, recursive = TRUE)
dir.create(unsafe_train, recursive = TRUE)
outside_file <- file.path(unsafe_root, "outside.txt")
writeLines("must survive", outside_file)
unsafe_audit <- file.path(unsafe_root, "cleanup.tsv")
writeLines("stale-success", unsafe_audit)
unsafe_owned <- artifact_rows(
  list("GCA_UNSAFE.1", "", "genome", "file", "../outside.txt")
)
unsafe_error <- tryCatch(
  {
    cleanup_owned_artifacts(
      unsafe_owned, data.frame(), unsafe_genomes, unsafe_train,
      unsafe_audit
    )
    NULL
  },
  error = function(error) error
)
check(
  inherits(unsafe_error, "error") && file.exists(outside_file),
  "unsafe relative paths fail without deleting outside the configured root"
)
unsafe_log <- read.delim(unsafe_audit, sep = "\t", stringsAsFactors = FALSE)
check(
  identical(
    names(unsafe_log),
    c("accessions", "artifact_type", "path", "action", "reason")
  ) && nrow(unsafe_log) == 1L && unsafe_log$action == "failed" &&
    grepl("Unsafe cleanup path", unsafe_log$reason, fixed = TRUE),
  "preflight failures atomically replace stale cleanup audits with the current error"
)

selected_unsafe_audit <- file.path(unsafe_root, "selected-cleanup.tsv")
selected_unsafe_error <- tryCatch(
  {
    cleanup_owned_artifacts(
      unsafe_owned,
      data.frame(
        out_acc = "GCA_UNSAFE.1", in1_acc = "GCA_SELECTED_IN1.1",
        in2_acc = "GCA_SELECTED_IN2.1", stringsAsFactors = FALSE
      ),
      unsafe_genomes, unsafe_train, selected_unsafe_audit
    )
    NULL
  },
  error = function(error) error
)
selected_unsafe_log <- read.delim(
  selected_unsafe_audit, sep = "\t", stringsAsFactors = FALSE
)
check(
  inherits(selected_unsafe_error, "error") &&
    selected_unsafe_log$action == "failed" &&
    grepl("Unsafe cleanup path", selected_unsafe_log$reason, fixed = TRUE),
  "a selected artifact that fails preflight is recorded as failed, not retained"
)

symlink_root <- tempfile("cleanup-symlink-")
symlink_genomes <- file.path(symlink_root, "genomes")
symlink_train <- file.path(symlink_root, "train")
dir.create(symlink_genomes, recursive = TRUE)
dir.create(symlink_train, recursive = TRUE)
symlink_target <- file.path(symlink_root, "target.txt")
symlink_path <- file.path(symlink_genomes, "owned-link")
writeLines("target survives", symlink_target)
invisible(file.symlink(symlink_target, symlink_path))
symlink_owned <- artifact_rows(
  list("GCA_LINK.1", "", "genome", "symlink", "owned-link")
)
cleanup_owned_artifacts(
  symlink_owned, data.frame(), symlink_genomes, symlink_train,
  file.path(symlink_root, "cleanup.tsv")
)
remaining_link <- Sys.readlink(symlink_path)
check(
  (is.na(remaining_link) || !nzchar(remaining_link)) && file.exists(symlink_target),
  "owned symlinks are unlinked without following their targets"
)

failure_root <- tempfile("cleanup-failure-")
failure_genomes <- file.path(failure_root, "genomes")
failure_train <- file.path(failure_root, "train")
dir.create(file.path(failure_genomes, "nonempty"), recursive = TRUE)
dir.create(failure_train, recursive = TRUE)
writeLines("pre-existing child", file.path(failure_genomes, "nonempty", "child.txt"))
failure_owned <- artifact_rows(
  list("GCA_FAILURE.1", "", "genome", "directory", "nonempty")
)
cleanup_error <- tryCatch(
  {
    cleanup_owned_artifacts(
      failure_owned, data.frame(), failure_genomes, failure_train,
      file.path(failure_root, "cleanup.tsv")
    )
    NULL
  },
  error = function(error) error
)
check(
  inherits(cleanup_error, "error") && dir.exists(file.path(failure_genomes, "nonempty")),
  "incomplete deletion is fatal and preserves an unowned child"
)

if (fail > 0L) {
  cat("\n", fail, " test(s) FAILED\n", sep = "")
  quit(status = 1)
}

cat("\nAll trio-selection helper tests passed.\n")
