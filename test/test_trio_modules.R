#!/usr/bin/env Rscript
# Loading and search-boundary regressions for the modular selector.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))
entry <- file.path(root, "src", "select", "trio_selection.R")
old_wd <- getwd()
setwd(tempdir())
source(entry)
stopifnot(identical(SCRIPT_DIR, dirname(entry)), exists("select_trios"))
e <- new.env(parent = globalenv())
e$run <- function(...) stop("external work during source")
source(entry, local = e)
setwd(old_wd)
stopifnot(identical(environment(e$select_trios), e),
          identical(environment(e$rank_assembly_candidates), e),
          identical(environment(e$record_selection_execution), e))
e$run <- function(...) "mock command"
stopifnot(identical(eval(quote(run("unused", character())), e), "mock command"))
# A moved function must resolve overrides in its source environment.
e$write_tsv <- function(...) NULL
e$read_tsv <- function(...) data.frame(excluded = "FALSE")
stopifnot(identical(e$evaluate_one(data.frame(), e$defaults, "unused", "unused")$excluded, "FALSE"))

e$candidate_outgroups <- function(...) c("skip", "missing", "pass", "later")
leaves <- data.frame(tip = c("a", "b", "skip", "missing", "pass", "later"),
                     accession = c("a", "b", "skip", "missing", "pass", "later"),
                     genus = c("A", "B", "A", "C", "D", "E"))
e$evaluate_one <- function(row, ...) {
  row$excluded <- "FALSE"
  row
}
e$build_row <- function(infoO, infoA, infoB, idt_12, idt_13, idt_23, tries) {
  data.frame(out_acc = infoO$accession, outgroup_tries = tries)
}
calls <- character()
identity <- function(a, b) {
  calls <<- c(calls, paste(a, b))
  if (a == "missing") NA_real_ else 95
}
opts <- e$defaults
opts$ingroup_idt_threshold <- 80
opts$max_outgroup_tries <- 2L
result <- e$select_trio_for_pair(NULL, leaves, c("a", "b"), opts, identity, "unused", "unused")
stopifnot(result$status == "selected", result$chosen$out_acc == "pass",
          result$tries == 2L, length(calls) == 5L,
          !any(grepl("skip|later", calls)))
opts$max_outgroup_tries <- 1L
calls <- character()
result <- e$select_trio_for_pair(NULL, leaves, c("a", "b"), opts, identity, "unused", "unused")
stopifnot(result$status == "no_valid_outgroup", result$tries == 1L, length(calls) == 3L)
cat("All trio module loading and search-boundary tests passed.\n")

e$write_tsv <- write_tsv
e$read_tsv <- read_tsv

# Partial cleanup failure must preserve completed deletions in the audit.
fixture <- tempfile("trio-partial-cleanup-")
dir.create(fixture)
paths <- file.path(fixture, c("first", "blocked", "later", "retained"))
file.create(paths)
artifacts <- data.frame(accession_1 = c("first", "blocked", "later", "keep"),
                        accession_2 = "", storage = "genome", artifact_type = "file",
                        relative_path = basename(paths))
selected <- data.frame(out_acc = "keep", in1_acc = "keep", in2_acc = "keep")
e$unlink <- function(x, ...) {
  if (identical(x, paths[2])) return(1L)
  base::unlink(x, ...)
}
audit_path <- file.path(fixture, "audit.tsv")
error <- try(e$cleanup_owned_artifacts(artifacts, selected, fixture, fixture, audit_path), silent = TRUE)
audit <- read.delim(audit_path)
stopifnot(inherits(error, "try-error"),
          identical(audit$action, c("deleted", "failed", "failed", "retained")),
          !file.exists(paths[1]), all(file.exists(paths[2:4])))
base::unlink(fixture, recursive = TRUE)
cat("Partial cleanup audit regression passed.\n")
