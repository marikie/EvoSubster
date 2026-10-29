# Helpers loaded by trio_selection.R into its calling environment.

IDENTITY_PREFIX <- "# substitution percent identity:"

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

remove_invalid_train_cache <- function(train) {
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
}

trained_owned_artifacts <- function(train, db_path, db_relative, db_preexisting, a1, a2) {
  records <- empty_owned_artifacts()
  train_link <- Sys.readlink(train)
  train_is_link <- !is.na(train_link) && nzchar(train_link)
  if (file.exists(train) && !dir.exists(train) && !train_is_link) {
    records <- rbind(records, data.frame(
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
    records <- rbind(records, data.frame(
      accession_1 = a1,
      accession_2 = "",
      storage = "train_cache",
      artifact_type = "directory_tree",
      relative_path = db_relative,
      stringsAsFactors = FALSE
    ))
  }
  records
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

  cache_dir <- resolve_train_cache_dir(opts)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(THREAD_NUM = opts$threads)
  run("python3", c(shQuote(file.path(SRC_DIR, "lib", "last_options.py")),
                   "guard", shQuote(cache_dir), "training"))

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
      remove_invalid_train_cache(train)
      fa1 <- get_fasta(a1); fa2 <- get_fasta(a2)
      db_relative <- paste0(a1, "db_", opts$date)
      db_path <- file.path(cache_dir, db_relative)
      db_preexisting <- path_exists_without_following(db_path)
      run("bash", c(shQuote(LAST_TRAIN), opts$date, shQuote(fa1), shQuote(fa2),
                    a1, a2, shQuote(cache_dir)))
      counters$trains <- counters$trains + 1L
      val <- parse_train_identity(train)
      if (!opts$keep_unused_species_data) {
        append_owned(trained_owned_artifacts(train, db_path, db_relative,
                                             db_preexisting, a1, a2))
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
