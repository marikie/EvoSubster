shape_manifest <- function(command, path, helper) {
  value <- system2("python3", c(shQuote(helper), command, shQuote(path)), stdout = TRUE)
  if (!is.null(attr(value, "status")) && attr(value, "status") != 0) stop("Manifest validation failed: ", path)
  jsonlite::fromJSON(paste(value, collapse = "\n"), simplifyVector = FALSE)
}
shape_ingest <- function(input, pattern, levels, helper, method) {
  runs <- unlist(shape_manifest("discover", input, helper), use.names = FALSE)
  vectors <- list()
  records <- list()
  for (run in sort(runs)) {
    m <- shape_manifest("read", run, helper)
    for (a in m$artifacts) {
      if (a$kind != if (pattern == "*_ncds.tsv") "single_ncds" else "single") next
      entries <- Filter(function(e) identical(e$slot, a$slot) && e$role != "outgroup", m$organisms)
      if (length(entries) != 1) next
      e <- entries[[1]]
      id <- paste(m$trio_key, m$run_key, a$slot, sep = "::")
      if (id %in% names(records)) stop("Duplicate sample identity: ", id)
      species <- e$raw_organism_name
      if (is.null(species)) species <- e$ncbi_full_name
      if (is.null(species)) species <- e$short_name
      rec <- list(
        sample_id = id, trio = m$trio_key, run = m$run_key, slot = a$slot, species = species,
        accession = if (is.null(e$accession)) "local" else e$accession,
        source = normalizePath(file.path(run, a$path)), included = FALSE, reason = ""
      )
      tax <- list()
      if (!is.null(e$taxonomy_json) && file.exists(file.path(run, e$taxonomy_json))) {
        t <- jsonlite::fromJSON(file.path(run, e$taxonomy_json), simplifyVector = FALSE)
        if (length(t$reports)) tax <- t$reports[[1]]$taxonomy$classification
      }
      for (level in levels) rec[[if (level == "species") "taxonomy_species" else level]] <- if (is.null(tax[[level]]$name)) "Unknown" else tax[[level]]$name
      value <- tryCatch(
        {
          x <- shape_vector(read.delim(rec$source, check.names = FALSE))
          if (method == "correlation" && sd(x) == 0) stop("Correlation undefined for constant spectrum")
          x
        },
        error = function(err) {
          rec$reason <<- conditionMessage(err)
          NULL
        }
      )
      if (!is.null(value)) {
        vectors[[id]] <- value
        rec$included <- TRUE
      }
      records[[id]] <- as.data.frame(rec, stringsAsFactors = FALSE)
    }
  }
  records <- if (length(records)) {
    do.call(rbind, records)
  } else {
    x <- as.data.frame(setNames(replicate(10 + length(levels), character(), simplify = FALSE), c("sample_id", "trio", "run", "slot", "species", "accession", "source", "included", "reason", "display_id", ifelse(levels == "species", "taxonomy_species", levels))))
    x
  }
  records$included <- as.logical(records$included)
  if (nrow(records)) records$display_id <- sprintf("S%04d (%s)", seq_len(nrow(records)), records$slot)
  matrix <- if (length(vectors)) do.call(rbind, vectors) else matrix(numeric(), 0, 96, dimnames = list(NULL, shape_channels()))
  list(matrix = matrix, records = records, runs = runs)
}
shape_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
shape_matrix_csv <- function(x, path) shape_csv(data.frame(sample_id = rownames(x), x, check.names = FALSE), path)
