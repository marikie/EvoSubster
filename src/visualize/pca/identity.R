# Bridge identity evaluation to the shared Python filtering helpers.

parse_identity_percentage <- function(value) {
  sanitized <- gsub("%", " ", value)
  matches <- stringr::str_extract_all(sanitized, "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?")[[1]]
  if (length(matches) == 0) {
    return(NA_real_)
  }
  as.numeric(matches[1])
}

collect_identity_strings <- function(run_dir, slot_map) {
  manifest <- attr(slot_map,"manifest") %||% manifest_command("read", run_dir)
  identity <- list()
  for (a in Filter(function(a) identical(a$kind, "train"), manifest$artifacts)) {
    key <- paste0("idt_", sub("org([123])-org([123])", "\\1\\2", a$slot))
    identity[[key]] <- paste(a$identity_percent, "%")
  }
  list(identity=identity, issues=character())
}

evaluate_dataset_filter <- function(slot_map, run_dir, options) {
  manifest <- attr(slot_map,"manifest") %||% manifest_command("read", run_dir)
  identities <- list()
  for (artifact in manifest$artifacts) {
    if (identical(artifact$kind, "train")) {
      key <- paste0("idt_", sub("org([123])-org([123])", "\\1\\2", artifact$slot))
      identities[[key]] <- artifact$identity_percent
    }
  }
  helper <- file.path(dirname(.pca_manifest_helper), "pca_identity_filter.py")
  payload <- jsonlite::toJSON(
    list(options = options, identities = identities, organisms = slot_map),
    auto_unbox = TRUE, null = "null", digits = NA
  )
  result <- system2("python3", shQuote(helper), input=payload, stdout=TRUE)
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop("Identity evaluation failed")
  }
  record <- jsonlite::fromJSON(paste(result,collapse="\n"), simplifyVector=FALSE)
  record$thresholds <- NULL
  for (slot in names(record$genera)) {
    record[[paste0("genus_", slot)]] <- record$genera[[slot]] %||% NA_character_
  }
  record$genera <- NULL
  record
}
