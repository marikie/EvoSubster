shape_region <- function(region = NULL, pattern = NULL) {
  mapped <- if (is.null(pattern)) NULL else if (pattern == "*.tsv") "whole" else "noncoding"
  if (!is.null(region) && !is.null(mapped) && region != mapped) stop("Conflicting --region and --tsv-pattern")
  if (!is.null(region)) region else if (!is.null(mapped)) mapped else "noncoding"
}

shape_sha256 <- function(path) {
  if (!file.exists(path)) {
    return("")
  }
  value <- system2("sha256sum", shQuote(path), stdout = TRUE)
  sub(" .*", "", value[[1]])
}

shape_legacy_ingest <- function(input, region, levels, method, run_pattern = "*") {
  vectors <- list()
  records <- list()
  decisions <- list()
  sources <- list()
  trios <- list.dirs(input, recursive = FALSE, full.names = TRUE)
  trios <- sort(unique(normalizePath(trios, mustWork = TRUE)))
  add_record <- function(trio, run, entry = NULL, reason = "", candidates = character()) {
    slot <- if (is.null(entry$slot)) "unknown" else entry$slot
    id <- paste(trio, basename(run), slot, sep = "::")
    species <- entry$raw_organism_name
    if (is.null(species)) species <- entry$ncbi_full_name
    if (is.null(species)) species <- entry$short_name
    if (is.null(species)) species <- "Unknown"
    rec <- list(
      sample_id = id, trio = trio, run = basename(run), slot = slot,
      species = species, accession = if (is.null(entry$accession)) "local" else entry$accession,
      source = paste(candidates, collapse = ";"), included = FALSE, reason = reason,
      region = region, source_sha256 = "", metadata_source = "", metadata_sha256 = ""
    )
    manifest <- file.path(run, "metadata", "metadata_manifest.json")
    rec$metadata_source <- manifest
    rec$metadata_sha256 <- shape_sha256(manifest)
    tax_path <- if (!is.null(entry$short_name) && !is.null(entry$accession)) {
      file.path(run, "metadata", paste0("taxonomy_", entry$short_name, "_", entry$accession, ".json"))
    } else {
      ""
    }
    tax <- tryCatch(
      {
        if (!nzchar(tax_path) || !file.exists(tax_path)) {
          list()
        } else {
          jsonlite::fromJSON(tax_path, simplifyVector = FALSE)$reports[[1]]$taxonomy$classification
        }
      },
      error = function(e) list()
    )
    rec$taxonomy_source <- if (file.exists(tax_path)) tax_path else ""
    rec$taxonomy_sha256 <- if (nzchar(rec$taxonomy_source)) shape_sha256(tax_path) else ""
    for (level in levels) {
      rec[[if (level == "species") "taxonomy_species" else level]] <-
        if (is.null(tax[[level]]$name)) "Unknown" else tax[[level]]$name
    }
    if (!nzchar(reason)) {
      rec$source_sha256 <- shape_sha256(candidates[[1]])
      value <- tryCatch(
        {
          x <- shape_vector(read.delim(candidates[[1]], check.names = FALSE))
          if (method == "correlation" && sd(x) == 0) stop("Correlation undefined for constant spectrum")
          x
        },
        error = function(e) {
          rec$reason <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(value)) {
        vectors[[id]] <<- value
        rec$included <- TRUE
      }
    }
    records[[id]] <<- as.data.frame(rec, stringsAsFactors = FALSE)
  }
  for (trio_path in trios) {
    trio <- trio_path
    runs <- list.dirs(trio_path, recursive = FALSE, full.names = TRUE)
    runs <- runs[grepl(glob2rx(run_pattern), basename(runs))]
    if (!length(runs)) next
    labels <- basename(runs)
    dates <- sub(".*?([0-9]{8}).*", "\\1", labels, perl = TRUE)
    parsed <- suppressWarnings(as.Date(dates, "%Y%m%d"))
    valid <- grepl("[0-9]{8}", labels) & !is.na(parsed)
    reason <- ""
    if (!all(valid)) reason <- "Run date unavailable or invalid"
    latest <- if (all(valid)) runs[parsed == max(parsed)] else runs
    if (!nzchar(reason) && length(latest) != 1) reason <- "Ambiguous latest run date"
    selected <- if (!nzchar(reason)) latest[[1]] else ""
    decisions[[length(decisions) + 1L]] <- data.frame(
      trio = trio, candidates = paste(runs, collapse = ";"),
      selected_run = selected, reason = if (nzchar(reason)) reason else "Latest matching run date", stringsAsFactors = FALSE
    )
    if (nzchar(reason)) {
      add_record(trio, trio_path, reason = reason)
      next
    }
    sources <- c(sources, list(selected))
    manifest <- file.path(selected, "metadata", "metadata_manifest.json")
    m <- tryCatch(jsonlite::fromJSON(manifest, simplifyVector = FALSE), error = function(e) NULL)
    object <- function(value) is.list(value) && !is.null(names(value))
    text_or_null <- function(value) is.null(value) || (is.character(value) && length(value) == 1L && !is.na(value))
    valid_manifest <- object(m) && is.list(m$organisms) && length(m$organisms) > 0L
    if (valid_manifest) {
      valid_manifest <- all(vapply(m$organisms, function(e) {
        object(e) && all(vapply(c("role", "slot", "short_name", "accession", "raw_organism_name", "ncbi_full_name"),
          function(field) text_or_null(e[[field]]), logical(1)))
      }, logical(1)))
    }
    if (!valid_manifest) {
      add_record(trio, selected, reason = "Missing or invalid legacy metadata manifest")
      next
    }
    entries <- Filter(function(e) identical(e$role, "ingroup"), m$organisms)
    if (!length(entries)) {
      add_record(trio, selected, reason = "No identifiable ingroups")
      next
    }
    slots <- vapply(entries, function(e) if (is.null(e$slot)) "" else e$slot, character(1))
    if (any(!slots %in% c("org2", "org3")) || anyDuplicated(slots)) {
      add_record(trio, selected, reason = "Ambiguous or invalid ingroup slots")
      next
    }
    for (e in entries) {
      if (is.null(e$slot) || is.null(e$short_name) || is.null(e$accession)) {
        add_record(trio, selected, e, "Incomplete ingroup identity")
        next
      }
      files <- list.files(file.path(selected, "statistics", "singlenuc"), pattern = "\\.tsv$", full.names = TRUE)
      names <- basename(files)
      match <- startsWith(names, paste0(e$accession, "_", e$short_name, "_")) |
        startsWith(names, paste0("cds_from_", e$short_name, "_"))
      match <- match & !grepl("dinuc", names, fixed = TRUE) &
        (grepl("ncds", names, fixed = TRUE) == (region == "noncoding"))
      candidates <- sort(files[match])
      reason <- if (!length(candidates)) "No TSV for requested region" else if (length(candidates) > 1) "Ambiguous TSV candidates" else ""
      add_record(trio, selected, e, reason, candidates)
    }
  }
  empty <- !length(records)
  if (empty) add_record("", input, reason = "No matching legacy runs found")
  records <- do.call(rbind, records)
  if (empty) records <- records[FALSE, , drop = FALSE]
  records$display_id <- sprintf("S%04d (%s)", seq_len(nrow(records)), records$slot)
  x <- if (length(vectors)) do.call(rbind, vectors) else matrix(numeric(), 0, 96, dimnames = list(NULL, shape_channels()))
  list(matrix = x, records = records, runs = unlist(sources), run_decisions = if (length(decisions)) do.call(rbind, decisions) else data.frame(trio = character(), candidates = character(), selected_run = character(), reason = character()))
}
