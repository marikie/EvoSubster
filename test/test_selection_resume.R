#!/usr/bin/env Rscript
# Exercise standalone selector persistence without downloads, training or plotting.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))
e <- new.env(parent = globalenv())
source(file.path(root, "src/select/trio_selection.R"), local = e)
work <- tempfile("selection resume ")
dir.create(work)
e$load_runtime_dependencies <- function() NULL
e$read.tree <- function(...) NULL
e$Ntip <- function(...) 3L
e$prune_to_best_assembly <- function(...) list(tree = NULL, leaves = data.frame())
e$write_stage0_tree <- function(...) c(nwk = "unused", pdf = "unused")
e$remove_stage1_tree_outputs <- function(...) NULL
e$write_selected_trio_tree <- function(...) c(nwk = "unused", pdf = "unused")
e$make_fetchers <- function(...) list(get_owned_artifacts = function() data.frame())
e$select_trios <- function(...) data.frame(out_acc = "GCA_1.1", in1_acc = "GCA_2.1", in2_acc = "GCA_3.1")
e$cleanup_owned_artifacts <- function(...) NULL
e$enumerate_pairs_dry <- function(...) data.frame()
args <- c("--tree", "unused.nwk", "--date", "20260818", "--threads", "2",
          "--out-dir", work, "--genome-dir", file.path(work, "genomes"))
record <- file.path(work, "selection_run.json")
e$main(args)
data <- jsonlite::fromJSON(record)
rows <- e$run("python3", c(shQuote(file.path(root, "src/lib/selection_run.py")),
                         "trios", shQuote(file.path(work, "selected_trios.tsv"))))
stopifnot(identical(rows, "GCA_1.1\tGCA_2.1\tGCA_3.1"))
stopifnot(data$date == "20260818", data$threads == 2,
          data$out_dir == dirname(normalizePath(work)),
          data$train_cache_dir == file.path(normalizePath(work), "train_cache"))
e$cleanup_owned_artifacts <- function(...) stop("cleanup failed")
error <- tryCatch(e$main(args), error = identity)
stopifnot(inherits(error, "error"), !file.exists(record))
writeLines("old record", record)
e$main(c(args, "--dry-run"))
stopifnot(!file.exists(record))
e$cleanup_owned_artifacts <- function(...) NULL
e$select_trios <- function(...) data.frame()
e$main(args)
stopifnot(file.exists(record))
rows <- e$run("python3", c(shQuote(file.path(root, "src/lib/selection_run.py")),
                         "trios", shQuote(file.path(work, "selected_trios.tsv"))))
stopifnot(length(rows) == 0)
unlink(work, recursive = TRUE)
cat("ok: standalone selection records settings only after success, including empty selection\n")
