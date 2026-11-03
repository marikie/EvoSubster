# Helpers loaded by trio_selection.R into its calling environment.

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


resolve_train_cache_dir <- function(opts) {
  if (!is.null(opts$train_cache_dir) && nzchar(opts$train_cache_dir)) {
    opts$train_cache_dir
  } else {
    file.path(opts$out_dir, "train_cache")
  }
}
