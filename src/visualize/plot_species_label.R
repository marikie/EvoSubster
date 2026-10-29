# Resolve the focal species from the run metadata without requiring network access.
plot_species_name <- function(input_path, full_name = NULL) {
  clean_name <- function(value) {
    if (is.null(value) || length(value) != 1 || is.na(value)) return("")
    trimws(gsub("_", " ", value, fixed = TRUE))
  }
  explicit_name <- clean_name(full_name)
  if (nzchar(explicit_name)) return(explicit_name)

  input_path <- normalizePath(input_path, mustWork = TRUE)
  filename <- basename(input_path)
  accession <- regmatches(filename, regexpr("^GC[AF]_[0-9]+[.][0-9]+", filename))
  short_name <- if (length(accession) && nzchar(accession)) {
    sub("_.*$", "", substring(filename, nchar(accession) + 2))
  } else {
    basename(dirname(dirname(input_path)))
  }

  directory <- dirname(input_path)
  repeat {
    manifest_path <- file.path(directory, "metadata", "metadata_manifest.json")
    if (file.exists(manifest_path)) {
      manifest <- tryCatch(
        jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
        error = function(error) {
          warning(sprintf("Could not read species metadata: %s", conditionMessage(error)))
          NULL
        }
      )
      organisms <- manifest$organisms
      matches <- Filter(function(entry) identical(entry$short_name, short_name), organisms)
      if (length(accession) && nzchar(accession)) {
        accession_matches <- Filter(function(entry) identical(entry$accession, accession), organisms)
        matches <- Filter(function(entry) identical(entry$short_name, short_name), accession_matches)
        if (length(matches) == 0) matches <- accession_matches
      }
      if (length(matches) == 1) {
        entry <- matches[[1]]
        for (candidate in c(entry$raw_organism_name, entry$ncbi_full_name, entry$directory_name)) {
          name <- clean_name(candidate)
          if (nzchar(name)) return(name)
        }
      }
      break
    }
    parent <- dirname(directory)
    if (identical(parent, directory)) break
    directory <- parent
  }
  warning(sprintf("Species full name unavailable for %s; supply it as the third argument.", filename))
  ""
}

add_species_label <- function(full_name, cex = 0.85, family = "sans") {
  if (nzchar(full_name)) {
    mtext(full_name, side = 3, line = 0.5, adj = 1, cex = cex, family = family, font = 3)
  }
}
