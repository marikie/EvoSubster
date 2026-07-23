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
    return(paste("Genus_species", fake_fasta, "", "", "", "", sep = "|"))
  }
  if (grepl("last_train[.]sh", args[1])) {
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

if (fail > 0L) {
  cat("\n", fail, " test(s) FAILED\n", sep = "")
  quit(status = 1)
}

cat("\nAll trio-selection helper tests passed.\n")
