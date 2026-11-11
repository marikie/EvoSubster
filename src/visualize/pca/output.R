# Write PCA tables, loadings, and the persistent run audit.

write_pc_loadings <- function(rotation_matrix, feature_names, output_dir, prefix) {
  file_path <- file.path(output_dir, sprintf("pc_loadings_%s.txt", prefix))
  sections <- list()
  format_section <- function(title, values) {
    formatted <- paste0(names(values), ": ", sprintf("%.6f", values))
    c(title, formatted)
  }
  if ("PC1" %in% colnames(rotation_matrix)) {
    pc1_values <- rotation_matrix[, "PC1"]
    names(pc1_values) <- feature_names
    pc1_sorted <- sort(pc1_values, decreasing = TRUE)
    sections <- c(
      sections,
      list(format_section("PC1 loadings (top 20):", head(pc1_sorted, 20))),
      list(format_section("PC1 loadings (bottom 20):", tail(pc1_sorted, 20)))
    )
  }
  if ("PC2" %in% colnames(rotation_matrix)) {
    pc2_values <- rotation_matrix[, "PC2"]
    names(pc2_values) <- feature_names
    pc2_sorted <- sort(pc2_values, decreasing = TRUE)
    sections <- c(
      sections,
      list(format_section("PC2 loadings (top 20):", head(pc2_sorted, 20))),
      list(format_section("PC2 loadings (bottom 20):", tail(pc2_sorted, 20)))
    )
  }
  if (length(sections) > 0) {
    content <- paste(unlist(sections), collapse = "\n")
    writeLines(content, file_path)
  }
}

write_empty_results_csv <- function(csv_path, classification_levels) {
  cols <- c("file", "accession", "short_name", "date", "metadata_genus", "species", "dataset", "cluster", "PC1", "PC2")
  for (level in classification_levels) cols <- c(cols, sprintf("classification_%s", level))
  df <- as.data.frame(matrix(ncol = length(cols), nrow = 0))
  names(df) <- cols
  readr::write_csv(df, csv_path)
}

save_results <- function(csv_path, files, pca_data, cluster_labels, metadata_records, classification_levels) {
  rows <- purrr::imap(metadata_records, function(rec, idx) {
    row <- list(
      file = basename(files[idx]),
      accession = rec$accession,
      short_name = rec$short_name,
      date = rec$date,
      metadata_genus = rec$genus,
      species = rec$species,
      dataset = rec$dataset,
      cluster = as.integer(cluster_labels[idx]),
      PC1 = pca_data[idx, 1],
      PC2 = pca_data[idx, 2]
    )
    for (level in classification_levels) {
      value <- rec$classifications[[level]] %||% "Unknown"
      row[[sprintf("classification_%s", level)]] <- value
    }
    row
  })
  df <- dplyr::bind_rows(rows)
  readr::write_csv(df, csv_path)
}

new_pca_audit <- function(settings, output_dir) {
  audit <- new.env(parent = emptyenv())
  audit$settings <- settings
  audit$records <- list()
  audit$output_dir <- output_dir
  audit
}

save_pca_audit <- function(audit) {
  settings <- audit$settings
  records <- audit$records
  output_dir <- audit$output_dir
  settings$finished_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  settings$trios <- length(records)
  settings$included_trios <- sum(vapply(records, function(x) isTRUE(x$included), logical(1)))
  settings$usable_samples <- sum(vapply(records, function(x) x$usable_samples, integer(1)))
  jsonlite::write_json(
    settings, file.path(output_dir, "pca_run_settings.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA
  )
  columns <- c(
    "trio_id", "layout", "run_path", "idt_12", "idt_13", "idt_23", "gap_12", "gap_13",
    "genus_org1", "genus_org2", "genus_org3", "genus_condition", "gap_waived",
    "ingroup_idt_threshold", "outgroup_idt_threshold", "min_idt_gap_between_out_in",
    "included", "reason", "usable_samples", "identity_source"
  )
  rows <- lapply(records, function(record) {
    values <- setNames(lapply(columns, function(key) record[[key]] %||% NA), columns)
    as.data.frame(values, stringsAsFactors = FALSE)
  })
  table <- if (length(rows)) {
    dplyr::bind_rows(rows)
  } else {
    as.data.frame(setNames(rep(list(logical()), length(columns)), columns))
  }
  write.csv(table, file.path(output_dir, "pca_trio_filter.csv"), row.names = FALSE, na = "")
}
