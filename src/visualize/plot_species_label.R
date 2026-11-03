# Resolve the focal species from the run metadata without requiring network access.
plot_species_name <- function(input_path, full_name = NULL) {
  clean_name <- function(value) {
    if (is.null(value) || length(value) != 1 || is.na(value)) return("")
    trimws(gsub("_", " ", value, fixed = TRUE))
  }
  explicit_name <- clean_name(full_name)
  if (nzchar(explicit_name)) return(explicit_name)

  input_path <- normalizePath(input_path, mustWork = TRUE)

  directory <- dirname(input_path)
  repeat {
    manifest_path <- file.path(directory, "metadata", "metadata_manifest.json")
    if (file.exists(manifest_path)) {
      manifest <- jsonlite::fromJSON(manifest_path,simplifyVector=FALSE)
      if (!identical(manifest$schema_version,2L) && !identical(manifest$schema_version,2)) {
        stop("Schema-2 manifest required for plot species labels")
      }
      matches <- Filter(function(artifact) {
        identical(normalizePath(file.path(directory,artifact$path),mustWork=FALSE),input_path)
      },manifest$artifacts)
      if (length(matches)>1) stop("Ambiguous plot input artifact: ",input_path)
      if (length(matches)==1) {
        organisms <- Filter(function(entry) identical(entry$slot,matches[[1]]$slot),manifest$organisms)
        if (length(organisms)!=1) stop("Missing unique organism slot for plot input")
        entry <- organisms[[1]]
        for (candidate in c(entry$raw_organism_name,entry$ncbi_full_name)) {
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
  warning(sprintf("Species full name unavailable for %s; supply it as the third argument.", basename(input_path)))
  ""
}

add_species_label <- function(full_name, cex = 0.85, family = "sans") {
  if (nzchar(full_name)) {
    mtext(full_name, side = 3, line = 0.5, adj = 1, cex = cex, family = family, font = 3)
  }
}
