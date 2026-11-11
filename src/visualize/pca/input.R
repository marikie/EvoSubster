# Discover runs and collect validated samples in their original order.

manifest_command <- function(command, path) {
  output <- system2("python3", c(shQuote(.pca_manifest_helper), command, shQuote(path)), stdout=TRUE)
  if (!is.null(attr(output, "status")) && attr(output, "status") != 0) stop("Manifest validation failed: ", path)
  jsonlite::fromJSON(paste(output, collapse="\n"), simplifyVector=FALSE)
}

warn <- function(msg) message(sprintf("Warning: %s", msg))

discover_dataset_dirs <- function(input_root) unlist(manifest_command("discover", input_root), use.names=FALSE)
find_latest_run_dir <- function(dataset_dir) {
  runs <- unlist(manifest_command("discover",dataset_dir),use.names=FALSE)
  if (!length(runs)) return(NULL)
  runs[[1]]
}


resolve_path <- function(value, base_dir) {
  if (is.null(value) || is.na(value) || value == "") {
    return(NULL)
  }
  path <- path.expand(value)
  if (!grepl("^/", path)) path <- file.path(base_dir, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}


load_manifest_lookup <- function(manifest_path) {
  manifest <- fromJSON(manifest_path, simplifyVector = FALSE)
  run_dir <- dirname(dirname(manifest_path))
  manifest <- manifest_command("read", run_dir)

  lookup <- list()
  organisms <- manifest$organisms %||% list()
  for (entry in organisms) {
    role <- tolower(entry$role %||% "")
    if (role == "outgroup") next
    short_name <- entry$short_name
    accession <- entry$accession
    metadata_json <- entry$metadata_json
    if (is.null(short_name) || is.null(metadata_json)) {
      stop(sprintf("Incomplete organism entry in manifest %s", manifest_path))
    }
    metadata_path <- file.path(run_dir, entry$metadata_json)
    taxonomy_path <- if (!is.null(entry$taxonomy_json)) file.path(run_dir, entry$taxonomy_json) else NA_character_
    key <- paste(accession %||% "local", short_name, sep = "::")
    lookup[[key]] <- list(metadata_path = metadata_path, taxonomy_path = taxonomy_path, entry = entry)
  }

  if (length(lookup) == 0) stop(sprintf("No ingroup organisms found in manifest %s", manifest_path))

  slot_map <- list()
  for (entry in organisms) {
    slot <- entry$slot
    if (!is.null(slot)) slot_map[[slot]] <- entry
  }

  attr(slot_map,"manifest") <- manifest
  list(manifest = manifest, lookup = lookup, slot_map = slot_map)
}

filter_entry_substitution_tsvs <- function(paths, entry) {
  prefix <- paste0(entry$accession, "_", entry$short_name, "_")
  paths[vapply(basename(paths), function(file_name) {
    startsWith(file_name, prefix) && !grepl("_dinuc_", file_name, fixed = TRUE)
  }, logical(1))]
}

substitution_tsv_search_locations <- function(run_dir, pattern, slot_map) {
  ingroup_entries <- Filter(function(entry) {
    tolower(entry$role %||% "") != "outgroup" &&
      !is.null(entry$short_name) && entry$short_name != ""
  }, slot_map)
  nested_locations <- vapply(
    ingroup_entries,
    function(entry) file.path(run_dir, "statistics", entry$short_name, "singlenuc", pattern),
    character(1)
  )
  unique(c(unname(nested_locations), file.path(run_dir, pattern)))
}

collect_substitution_tsvs <- function(run_dir, pattern, slot_map) {
  if (!pattern %in% c("*_ncds.tsv","*.tsv")) stop("Unsupported --tsv-pattern; use *_ncds.tsv (noncoding) or *.tsv (all genomic).")
  kind <- if (grepl("ncds", pattern, fixed=TRUE)) "single_ncds" else "single"
  manifest <- attr(slot_map,"manifest") %||% manifest_command("read", run_dir)
  unlist(lapply(Filter(function(a) identical(a$kind, kind), manifest$artifacts),
                function(a) file.path(run_dir, a$path)), use.names=FALSE)
}

load_tsv_vectors <- function(tsv_path, reference_order = NULL) {
  df <- read_tsv(tsv_path, show_col_types = FALSE)
  if (is.null(reference_order)) {
    df_sorted <- df %>% arrange(.data$mutType)
    observed <- compute_observed_vector(df_sorted)
    list(
      logratio = compute_logratio_vector(df_sorted, observed),
      observed = observed,
      order = df_sorted$mutType
    )
  } else {
    df_indexed <- df %>%
      filter(.data$mutType %in% reference_order) %>%
      arrange(match(.data$mutType, reference_order))
    if (!all(reference_order %in% df_indexed$mutType)) {
      stop(sprintf("TSV file %s does not contain expected mutType values.", basename(tsv_path)))
    }
    observed <- compute_observed_vector(df_indexed)
    list(
      logratio = compute_logratio_vector(df_indexed, observed),
      observed = observed,
      order = reference_order
    )
  }
}

collect_pca_inputs <- function(input_root, args, options, audit) {
  metadata_cache <- list()
  taxonomy_cache <- list()
  logratio_vectors <- list()
  observed_vectors <- list()
  species_labels <- character()
  metadata_records <- list()
  files <- character()
  reference_order <- NULL
  observed_levels <- character()
  dataset_status <- list()

  dataset_dirs <- discover_dataset_dirs(input_root)

  for (dataset_dir in dataset_dirs) {
    dataset_name <- basename(dataset_dir)
    latest_run <- find_latest_run_dir(dataset_dir)
    if (is.null(latest_run)) {
      warn(sprintf("%s: No dated run directories found under %s.", dataset_name, dataset_dir))
      next
    }
    manifest_path <- file.path(latest_run, "metadata", "metadata_manifest.json")
    if (!file.exists(manifest_path)) {
      warn(sprintf("%s: metadata manifest not found at %s.", dataset_name, manifest_path))
      next
    }
    manifest_bundle <- load_manifest_lookup(manifest_path)
    if (is.null(manifest_bundle)) next
    dataset_name <- manifest_bundle$manifest$trio_key
    slot_map <- manifest_bundle$slot_map
    record <- evaluate_dataset_filter(slot_map, latest_run, options)
    record$trio_id <- dataset_name
    record$layout <- manifest_bundle$manifest$layout %||% "v2"
    record$run_path <- latest_run
    record$usable_samples <- 0L
    record$reason <- paste(Filter(nzchar, c(record$reason, unlist(manifest_bundle$manifest$input_issues))), collapse = "; ")
    audit$records[[dataset_name]] <- record
    dataset_status[[dataset_name]] <- record
    if (!isTRUE(record$included)) next

    add_sample_issue <- function(message) {
      audit$records[[dataset_name]]$reason <- paste(
        Filter(nzchar, c(audit$records[[dataset_name]]$reason, message)), collapse = "; "
      )
      warn(sprintf("%s: %s", dataset_name, message))
    }
    selected_kind <- if (identical(args$tsv_pattern, "*_ncds.tsv")) "single_ncds" else "single"
    selected_artifacts <- Filter(function(a) identical(a$kind, selected_kind), manifest_bundle$manifest$artifacts)
    for (slot in c("org2", "org3")) {
      if (!any(vapply(selected_artifacts, function(a) identical(a$slot, slot), logical(1)))) {
        add_sample_issue(sprintf("%s: Missing canonical %s TSV", slot, selected_kind))
      }
    }
    tsv_files <- collect_substitution_tsvs(latest_run, args$tsv_pattern, slot_map)
    if (length(tsv_files) == 0) {
      locations <- substitution_tsv_search_locations(latest_run, args$tsv_pattern, slot_map)
      warn(sprintf(
        "%s: No TSV files matching '%s' in %s. Checked: %s.",
        dataset_name,
        args$tsv_pattern,
        latest_run,
        paste(locations, collapse = "; ")
      ))
      next
    }
    samples_before <- length(logratio_vectors)
    for (tsv_file in sort(tsv_files)) {
      candidates <- Filter(function(a) identical(normalizePath(file.path(latest_run, a$path)), normalizePath(tsv_file)), manifest_bundle$manifest$artifacts)
      if (length(candidates) != 1) stop("Ambiguous TSV artifact: ", tsv_file)
      entry <- slot_map[[candidates[[1]]$slot]]
      accession <- entry$accession %||% "local"
      short_name <- entry$short_name
      date <- manifest_bundle$manifest$date
      key <- paste(entry$accession %||% "local", short_name, sep="::")
      manifest_entry <- manifest_bundle$lookup[[key]]
      metadata_path <- manifest_entry$metadata_path
      metadata_json <- tryCatch(load_json_cached(metadata_path, metadata_cache), error = function(e) {
        add_sample_issue(sprintf("Failed to load metadata JSON %s: %s", metadata_path, e$message))
        NULL
      })
      if (is.null(metadata_json)) next
      taxonomy_path <- manifest_entry$taxonomy_path
      taxonomy_json <- if (!is.na(taxonomy_path) && file.exists(taxonomy_path)) {
        tryCatch(load_json_cached(taxonomy_path, taxonomy_cache), error = function(e) {
          add_sample_issue(sprintf("Unusable taxonomy JSON %s: %s", taxonomy_path, e$message))
          list()
        })
      } else {
        add_sample_issue(sprintf("%s: Missing local taxonomy JSON", entry$slot))
        list()
      }
      genus_info <- if (is.null(entry$raw_organism_name) && is.null(entry$ncbi_full_name)) {
        list(genus="Unknown",species=short_name)
      } else extract_genus(metadata_json, short_name, accession)
      classification_map <- extract_classification_map(taxonomy_json, sprintf("%s (%s)", short_name, accession))
      observed_levels <- unique(c(observed_levels, names(classification_map)))

      vectors <- tryCatch(load_tsv_vectors(tsv_file, reference_order), error = function(e) {
        add_sample_issue(sprintf("Failed to parse TSV %s: %s", basename(tsv_file), e$message))
        NULL
      })
      if (is.null(vectors)) next
      reference_order <- vectors$order
      logratio_vectors[[length(logratio_vectors) + 1]] <- vectors$logratio
      observed_vectors[[length(observed_vectors) + 1]] <- vectors$observed
      species_labels <- c(species_labels, genus_info$species)
      files <- c(files, tsv_file)
      metadata_records[[length(metadata_records) + 1]] <- list(
        accession = accession,
        short_name = short_name,
        date = date,
        genus = genus_info$genus,
        species = genus_info$species,
        classifications = classification_map,
        dataset = dataset_name
      )
    }
    added <- length(logratio_vectors) - samples_before
    audit$records[[dataset_name]]$usable_samples <- as.integer(added)
    if (added == 0) warn(sprintf("%s: No usable TSV files after validation; dataset skipped.", dataset_name))
  }
  list(
    logratio_vectors = logratio_vectors,
    observed_vectors = observed_vectors,
    species_labels = species_labels,
    metadata_records = metadata_records,
    files = files,
    reference_order = reference_order,
    observed_levels = observed_levels,
    dataset_status = dataset_status
  )
}
