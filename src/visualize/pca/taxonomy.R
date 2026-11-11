# Read local metadata and derive taxonomy labels.

load_json_cached <- function(path, cache) {
  if (is.null(path) || !file.exists(path)) stop(sprintf("Required JSON file not found: %s", path))
  if (!path %in% names(cache)) cache[[path]] <- fromJSON(path, simplifyVector = FALSE)
  cache[[path]]
}

extract_tax_id <- function(metadata) {
  reports <- metadata$reports
  if (is.null(reports) || length(reports) == 0) {
    return(NULL)
  }
  first <- reports[[1]]
  tax_id <- first$organism$tax_id %||% first$biosample$description$organism$tax_id
  if (is.null(tax_id) || tax_id == "" || tax_id == "null") {
    return(NULL)
  }
  as.character(tax_id)
}

ensure_taxonomy_json <- function(taxonomy_path, metadata, short_name, accession) {
  if (!is.null(taxonomy_path) && file.exists(taxonomy_path) && file.size(taxonomy_path) > 0) {
    return(TRUE)
  }
  tax_id <- extract_tax_id(metadata)
  if (is.null(tax_id)) {
    warn(sprintf("%s: tax_id not found; cannot fetch taxonomy JSON.", short_name))
    return(FALSE)
  }
  dir.create(dirname(taxonomy_path), recursive = TRUE, showWarnings = FALSE)
  command <- sprintf('datasets summary taxonomy taxon "%s"', tax_id)
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  result <- tryCatch(system(command, intern = TRUE, ignore.stderr = TRUE), warning = function(w) "", error = function(e) "")
  if (length(result) == 0) {
    warn(sprintf("%s: Failed to download taxonomy for tax_id %s.", short_name, tax_id))
    return(FALSE)
  }
  writeLines(result, tmp)
  if (!file.rename(tmp, taxonomy_path)) {
    if (!file.copy(tmp, taxonomy_path, overwrite = TRUE)) {
      unlink(tmp)
      warn(sprintf("%s: Failed to write taxonomy JSON for tax_id %s.", short_name, tax_id))
      return(FALSE)
    }
    unlink(tmp)
  }
  TRUE
}

extract_classification_map <- function(taxonomy, identifier) {
  reports <- taxonomy$reports
  if (is.null(reports) || length(reports) == 0) {
    warn(sprintf("taxonomy JSON missing 'reports' for %s", identifier))
    return(list())
  }
  classification <- reports[[1]]$taxonomy$classification %||% list()
  labels <- list()
  for (level_name in names(classification)) {
    entry <- classification[[level_name]]
    if (!is.list(entry)) next
    name <- entry$name
    if (!is.null(name)) labels[[tolower(level_name)]] <- name
  }
  if (length(labels) == 0) warn(sprintf("taxonomy classification map missing for %s", identifier))
  labels
}

extract_genus <- function(metadata, short_name, accession) {
  reports <- metadata$reports
  if (!is.null(reports) && length(reports) > 0) {
    organism <- reports[[1]]$organism
    organism_name <- organism$organism_name
    if (!is.null(organism_name)) {
      genus <- strsplit(organism_name, " +")[[1]][1]
      return(list(genus = genus, species = organism_name))
    }
  }
  list(genus = "Unknown", species = accession)
}

extract_genus_label <- function(entry) {
  candidates <- c(entry$raw_organism_name, entry$ncbi_full_name)
  for (candidate in candidates) {
    if (is.null(candidate) || candidate == "") next
    sanitized <- trimws(gsub("_", " ", candidate))
    if (sanitized == "") next
    token <- strsplit(sanitized, " +")[[1]][1]
    if (!is.null(token) && token != "") {
      return(tolower(token))
    }
  }
  NULL
}

evaluate_genus_condition <- function(slot_map) {
  slots <- c("org1", "org2", "org3")
  genus_values <- list()
  for (slot in slots) {
    entry <- slot_map[[slot]]
    if (is.null(entry)) {
      return(list(result = NULL, issues = sprintf("%s: Missing metadata entry for genus evaluation.", slot)))
    }
    genus <- extract_genus_label(entry)
    if (is.null(genus)) {
      return(list(result = NULL, issues = sprintf("%s: Unable to derive genus from metadata.", slot)))
    }
    genus_values[[slot]] <- genus
  }
  condition <- genus_values$org1 != genus_values$org2 &&
    genus_values$org1 != genus_values$org3 &&
    genus_values$org2 == genus_values$org3
  list(result = condition, issues = character())
}
