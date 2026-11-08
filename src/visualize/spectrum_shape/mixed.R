shape_mixed_ingest <- function(input, region, levels, method, run_pattern, helper) {
  trios <- sort(unique(normalizePath(list.dirs(input, recursive = FALSE), mustWork = TRUE)))
  bundles <- list()
  decisions <- list()
  for (trio in trios) {
    runs <- list.dirs(trio, recursive = FALSE)
    runs <- runs[grepl(glob2rx(run_pattern), basename(runs))]
    if (!length(runs)) next
    times <- shape_run_time(basename(runs))
    if (anyNA(times) || sum(times == max(times)) != 1) stop("Ambiguous or invalid matching run dates: ", trio)
    alias <- runs[which.max(times)]
    real <- normalizePath(alias, mustWork = TRUE)
    manifest_path <- file.path(real, "metadata", "metadata_manifest.json")
    m <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
    v2 <- any(c("schema_version", "run_key", "trio_key", "artifacts", "status") %in% names(m))
    if (v2) {
      m <- shape_manifest("read", real, helper)
      if (!identical(m$status, "complete")) stop("Selected v2 run is not complete: ", real)
      # Compare ordered accession identities against older cohort metadata.
      older <- setdiff(runs, alias)
      older <- older[order(shape_run_time(basename(older)), decreasing = TRUE)]
      if (length(older)) {
        previous <- jsonlite::fromJSON(file.path(older[1], "metadata", "metadata_manifest.json"), simplifyVector = FALSE)
        accessions <- function(value) vapply(c("org1", "org2", "org3"), function(slot) {
          entries <- Filter(function(e) identical(e$slot, slot), value$organisms)
          if (length(entries) != 1) stop("Ambiguous organism identity")
          entries[[1]]$accession
        }, character(1))
        if (!identical(accessions(m), accessions(previous))) stop("Selected v2 accessions differ from cohort trio: ", trio)
      }
      bundle <- shape_ingest(real, if (region == "whole") "*.tsv" else "*_ncds.tsv", levels, helper, method)
      for (entry in Filter(function(e) identical(e$role, "ingroup"), m$organisms)) {
        if (entry$slot %in% bundle$records$slot) next
        row <- list(sample_id = paste(m$trio_key, m$run_key, entry$slot, sep = "::"),
          trio = m$trio_key, run = m$run_key, slot = entry$slot,
          species = if (is.null(entry$raw_organism_name)) entry$short_name else entry$raw_organism_name,
          accession = entry$accession, source = "", included = FALSE,
          reason = "No TSV for requested region", display_id = "")
        tax <- if (is.null(entry$taxonomy_json)) list() else jsonlite::fromJSON(file.path(real, entry$taxonomy_json), simplifyVector = FALSE)$reports[[1]]$taxonomy$classification
        for (level in levels) row[[if (level == "species") "taxonomy_species" else level]] <- if (is.null(tax[[level]]$name)) "Unknown" else tax[[level]]$name
        bundle$records <- rbind(bundle$records, as.data.frame(row, stringsAsFactors = FALSE))
      }
    } else {
      bundle <- shape_legacy_ingest(input, region, levels, method, run_pattern, setNames(list(alias), trio))
    }
    bundle$records$logical_trio <- trio
    bundle$records$selection_alias <- alias
    bundle$records$resolved_run <- real
    bundle$records$input_format <- if (v2) "v2" else "legacy"
    bundle$records$region <- region
    bundles[[length(bundles) + 1L]] <- bundle
    decisions[[length(decisions) + 1L]] <- data.frame(trio = trio, selected_run = alias, resolved_run = real,
      format = if (v2) "v2" else "legacy", reason = "Latest matching run timestamp")
  }
  if (!length(bundles)) stop("No matching runs found")
  columns <- unique(unlist(lapply(bundles, function(b) names(b$records))))
  records <- do.call(rbind, lapply(bundles, function(b) {
    for (column in setdiff(columns, names(b$records))) b$records[[column]] <- ""
    b$records[, columns, drop = FALSE]
  }))
  records$display_id <- sprintf("S%04d (%s)", seq_len(nrow(records)), records$slot)
  x <- do.call(rbind, lapply(bundles, function(b) b$matrix))
  list(matrix = x, records = records, runs = vapply(bundles, function(b) b$runs[[1]], character(1)),
       run_decisions = do.call(rbind, decisions))
}
