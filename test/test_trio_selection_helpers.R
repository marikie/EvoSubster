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
  "Mus_musculus_GCF_000001635.27",
  "Drosophila_melanogaster_strain_A_GCA_000001215.4"
)
check(
  identical(
    unname(species_from_label(labels)),
    c("Homo sapiens", "Mus musculus", "Drosophila melanogaster")
  ),
  "species parser accepts accession-less and accession-suffixed genus_species labels"
)

malformed <- species_from_label(c("Homo", "_", "GCA_000001405.40"))
check(
  all(is.na(malformed)),
  "species parser rejects labels without both genus and species tokens"
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
