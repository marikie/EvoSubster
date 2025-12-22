#!/usr/bin/env Rscript

# Walk dataset directories, gather substitution TSV files from each latest run,
# and run 2D PCA clustering (logratio + observed variants, normal + filtered).

suppressPackageStartupMessages({
  library(argparse)
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(cluster)
  library(scales)
  library(ggpubr)
  library(rlang)
  library(ggrepel)
})

EPS <- 1e-12
IDENTITY_PREFIX <- "# substitution percent identity:"
DEFAULT_TSV_PATTERN <- "*_maflinked_ncds.tsv"
CLASSIFICATION_ORDER <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

invisible(utils::globalVariables(c("mutNum", "totalRootNum", "mutType")))

parse_args <- function() {
  parser <- ArgumentParser(
    description = "Traverse dataset directories, collect TSV files, and run PCA clustering."
  )
  parser$add_argument("input_dir", help = "Root directory containing dataset subdirectories (e.g., data/cnidaria).")
  parser$add_argument("--tsv-pattern", default = DEFAULT_TSV_PATTERN, help = "Glob pattern for TSV files inside each run dir.")
  parser$add_argument("--output-dir", help = "Directory for PCA outputs (default: <input_dir>/pca).")
  parser$add_argument("--max-clusters", type = "integer", default = 10L, help = "Maximum clusters to evaluate.")
  parser$add_argument("--random-state", type = "integer", default = 42L, help = "Random seed.")
  parser$parse_args()
}

warn <- function(msg) message(sprintf("Warning: %s", msg))

find_latest_run_dir <- function(dataset_dir) {
  entries <- list.dirs(dataset_dir, full.names = TRUE, recursive = FALSE)
  runs <- entries[grepl("^[0-9]+$", basename(entries))]
  if (length(runs) == 0) {
    return(NULL)
  }
  runs[which.max(basename(runs))]
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
  metadata_dir_value <- manifest$metadata_dir
  if (is.null(metadata_dir_value)) stop(sprintf("'metadata_dir' missing in manifest %s", manifest_path))
  metadata_dir <- resolve_path(metadata_dir_value, dirname(manifest_path))
  if (is.null(metadata_dir) || !dir.exists(metadata_dir)) {
    stop(sprintf("metadata_dir does not exist: %s", metadata_dir_value))
}

  lookup <- list()
  organisms <- manifest$organisms %||% list()
  for (entry in organisms) {
    role <- tolower(entry$role %||% "")
    if (role == "outgroup") next
    short_name <- entry$short_name
    accession <- entry$accession
    metadata_json <- entry$metadata_json
    if (is.null(short_name) || is.null(accession) || is.null(metadata_json)) {
      stop(sprintf("Incomplete organism entry in manifest %s", manifest_path))
    }
    metadata_path <- resolve_path(metadata_json, dirname(manifest_path))
    taxonomy_path <- file.path(metadata_dir, sprintf("taxonomy_%s_%s.json", short_name, accession))
    key <- paste(accession, short_name, sep = "::")
    lookup[[key]] <- list(metadata_path = metadata_path, taxonomy_path = taxonomy_path, entry = entry)
  }

  if (length(lookup) == 0) stop(sprintf("No ingroup organisms found in manifest %s", manifest_path))

  slot_map <- list()
  for (entry in organisms) {
    slot <- entry$slot
    if (!is.null(slot)) slot_map[[slot]] <- entry
  }

  list(manifest = manifest, lookup = lookup, slot_map = slot_map)
}

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
  candidates <- c(entry$raw_organism_name, entry$ncbi_full_name, entry$directory_name, entry$short_name)
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

select_first <- function(run_dir, pattern) {
  matches <- Sys.glob(file.path(run_dir, pattern))
  if (length(matches) == 0) {
    return(NULL)
  }
  sort(matches)[1]
}

extract_identity_line <- function(train_path) {
  lines <- readLines(train_path, warn = FALSE)
  hits <- lines[startsWith(trimws(lines), IDENTITY_PREFIX)]
  if (length(hits) >= 2) {
    return(trimws(hits[length(hits) - 1]))
  }
  if (length(hits) == 1) {
    return(trimws(hits[1]))
  }
  NULL
}

extract_identity_value <- function(line) {
  parts <- strsplit(line, ":", fixed = TRUE)[[1]]
  trimws(tail(parts, 1))
}

parse_identity_percentage <- function(value) {
  sanitized <- gsub("%", " ", value)
  matches <- stringr::str_extract_all(sanitized, "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?")[[1]]
  if (length(matches) == 0) {
    return(NA_real_)
  }
  as.numeric(matches[1])
}

collect_identity_strings <- function(run_dir, slot_map) {
  identity <- list()
  issues <- character()
  short_name_for <- function(slot) slot_map[[slot]]$short_name %||% NA
  pair_specs <- list(
    list(key = "idt_12", a = "org1", b = "org2"),
    list(key = "idt_13", a = "org1", b = "org3"),
    list(key = "idt_23", a = "org2", b = "org3")
  )
  for (spec in pair_specs) {
    short_a <- short_name_for(spec$a)
    short_b <- short_name_for(spec$b)
    if (is.na(short_a) || is.na(short_b)) {
      issues <- c(issues, sprintf("%s: Missing short_name for %s or %s.", spec$key, spec$a, spec$b))
      next
    }
    pattern <- sprintf("*%s2%s_*.train", short_a, short_b)
    train_path <- select_first(run_dir, pattern)
    if (is.null(train_path)) {
      issues <- c(issues, sprintf("%s: No train file matching '%s'.", spec$key, pattern))
      next
    }
    identity_line <- extract_identity_line(train_path)
    if (is.null(identity_line)) {
      issues <- c(issues, sprintf("%s: '%s' line missing in %s.", spec$key, IDENTITY_PREFIX, basename(train_path)))
      next
    }
    raw_value <- extract_identity_value(identity_line)
    value <- if (endsWith(raw_value, "%")) raw_value else paste(raw_value, "%")
    identity[[spec$key]] <- value
  }
  list(identity = identity, issues = issues)
}

evaluate_identity_condition <- function(identity_values) {
  if (length(identity_values) == 0) {
    return(list(result = NULL, issues = "Identity values missing; cannot evaluate condition."))
  }
  numeric_values <- list()
  for (key in c("idt_12", "idt_13", "idt_23")) {
    raw_value <- identity_values[[key]]
    if (is.null(raw_value)) {
      return(list(result = NULL, issues = sprintf("%s: Identity value missing.", key)))
    }
    parsed <- parse_identity_percentage(raw_value)
    if (is.na(parsed)) {
      return(list(result = NULL, issues = sprintf("%s: Unable to parse numeric value from '%s'.", key, raw_value)))
    }
    numeric_values[[key]] <- parsed
  }
  condition <- numeric_values$idt_12 < numeric_values$idt_23 && numeric_values$idt_13 < numeric_values$idt_23
  list(result = condition, issues = character())
}

evaluate_dataset_filter <- function(slot_map, run_dir) {
  genus_res <- evaluate_genus_condition(slot_map)
  identity_values <- list()
  identity_issues <- character()
  if (!is.null(run_dir) && dir.exists(run_dir)) {
    identity_res <- collect_identity_strings(run_dir, slot_map)
    identity_values <- identity_res$identity
    identity_issues <- identity_res$issues
  } else {
    identity_issues <- "Identity evaluation skipped: run_dir missing or does not exist."
  }
  identity_condition <- list(result = NULL, issues = character())
  if (length(identity_values) > 0) {
    identity_condition <- evaluate_identity_condition(identity_values)
  }
  excluded <- identical(genus_res$result, FALSE) && identical(identity_condition$result, FALSE)
  list(
    filter_status = list(
      genus_condition_met = genus_res$result,
      identity_condition_met = identity_condition$result,
      excluded = excluded
    ),
    issues = c(genus_res$issues, identity_issues, identity_condition$issues)
  )
}

compute_observed_vector <- function(df) {
  if (!all(c("mutType", "mutNum", "totalRootNum") %in% names(df))) {
    stop("TSV file is missing required columns (mutType, mutNum, totalRootNum)")
  }
  totals <- as.numeric(df$totalRootNum)
  mut_counts <- as.numeric(df$mutNum)
  observed <- ifelse(totals != 0, mut_counts / totals, 0)
  observed
}

compute_logratio_vector <- function(df, observed = NULL) {
  if (is.null(observed)) observed <- compute_observed_vector(df)
  mean_observed <- ifelse(length(observed) > 0, mean(observed), 0)
  adjusted_mean <- max(mean_observed, EPS)
  ratios <- (observed + EPS) / adjusted_mean
  log2(ratios)
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

determine_clusters <- function(data, max_clusters, seed) {
  sample_count <- nrow(data)
  if (sample_count < 2) stop("At least two samples are required for clustering.")
  if (sample_count == 2) {
    return(list(best_k = 2, inertia = numeric(), silhouette = numeric()))
  }
  adjusted_max <- min(max_clusters, sample_count - 1)
  if (adjusted_max < 2) adjusted_max <- 2
  set.seed(seed)
  inertia <- numeric()
  silhouette_scores <- numeric()
  for (k in 2:adjusted_max) {
    model <- kmeans(data, centers = k, nstart = 10)
    inertia <- c(inertia, model$tot.withinss)
    silhouette_scores <- c(silhouette_scores, mean(cluster::silhouette(model$cluster, stats::dist(data))[, "sil_width"]))
  }
  best_k <- which.max(silhouette_scores) + 2 - 1
  list(best_k = best_k, inertia = inertia, silhouette = silhouette_scores)
}

plot_metrics <- function(metrics_path, inertia, silhouette) {
  if (length(inertia) == 0 && length(silhouette) == 0) {
    ggsave(metrics_path, ggplot() +
      theme_void(), width = 8, height = 4, dpi = 300)
    return()
  }
  p1 <- if (length(inertia) > 0) {
    df <- tibble(k = seq_along(inertia) + 1, inertia = inertia)
    ggplot(df, aes(x = .data$k, y = .data$inertia)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    labs(title = "Elbow Method", x = "Number of Clusters", y = "Inertia") +
      theme_minimal()
  } else {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient samples\nfor inertia plot", size = 4) +
      theme_void()
  }
  p2 <- if (length(silhouette) > 0) {
    df <- tibble(k = seq_along(silhouette) + 1, silhouette = silhouette)
    ggplot(df, aes(x = .data$k, y = .data$silhouette)) +
    geom_line(color = "red") +
    geom_point(color = "red") +
    labs(title = "Silhouette Scores", x = "Number of Clusters", y = "Silhouette Score") +
      theme_minimal()
  } else {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Silhouette undefined\n(n_samples too small)", size = 4) +
      theme_void()
  }
  ggsave(metrics_path, ggpubr::ggarrange(p1, p2, nrow = 1), dpi = 300, width = 12, height = 5)
}

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

plot_pca_scatter <- function(output_path, pca_data, species_labels, color_labels, cluster_labels, label_name, feature_desc) {
  df <- tibble(
    PC1 = pca_data[, 1],
    PC2 = pca_data[, 2],
    species = stringr::str_replace_all(species_labels, "_", " "),
    label = color_labels,
    cluster = factor(cluster_labels)
  )
  palette <- scales::hue_pal()(length(unique(df$label)))
  names(palette) <- unique(df$label)
  classification_label <- stringr::str_to_title(label_name)
  p <- ggplot(df, aes(x = .data$PC1, y = .data$PC2, color = .data$label)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_text_repel(
      aes(label = .data$species),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.35,
      max.overlaps = Inf,
      min.segment.length = 0,
      force = 4,
      force_pull = 0.1,
      max.iter = 5000,
      max.time = 5,
      seed = 42,
      segment.color = "#8c8c8c",
      segment.size = 0.25,
      fontface = "italic"
    ) +
    scale_color_manual(values = palette, name = classification_label) +
    labs(
      title = sprintf("K-means Clustering (k=%d)", length(unique(cluster_labels))),
      subtitle = sprintf("Features: %s", feature_desc),
      x = "PC1",
      y = "PC2"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "#dcdcdc", linewidth = 0.45),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_text(face = "plain"),
      legend.position = "right"
    )
  ggsave(output_path, p, dpi = 300, width = 10, height = 8)
}

write_placeholder_plot <- function(output_path, message_text) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = message_text, size = 4, hjust = 0.5, vjust = 0.5) +
    theme_void()
  ggsave(output_path, p, dpi = 300, width = 6, height = 4)
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

run_pca_for_indices <- function(vectors, feature_order, indices, output_prefix, plot_levels, metadata_records, files, classification_levels, species_labels, max_clusters, random_state, output_dir) {
  if (length(indices) < 2) stop(sprintf("%s PCA requires at least two samples.", output_prefix))
  subset_vectors <- vectors[indices]
  matrix_data <- do.call(rbind, subset_vectors)
  scaled_data <- scale(matrix_data)
  pca <- prcomp(scaled_data, center = FALSE, scale. = FALSE)
  pca_data <- pca$x[, 1:2, drop = FALSE]
  if (is.null(colnames(pca$rotation))) {
    colnames(pca$rotation) <- feature_names
  }
  write_pc_loadings(pca$rotation, feature_order, output_dir, output_prefix)
  cluster_info <- determine_clusters(pca_data, max_clusters, random_state)
  plot_metrics(file.path(output_dir, sprintf("clustering_metrics_%s.png", output_prefix)), cluster_info$inertia, cluster_info$silhouette)
  set.seed(random_state)
  sample_count <- nrow(pca_data)
  message(sprintf("Running %s PCA: samples=%d, best_k=%s", output_prefix, sample_count, cluster_info$best_k))
  if (cluster_info$best_k >= sample_count) {
    cluster_labels <- seq_len(sample_count)
  } else {
  final_model <- kmeans(pca_data, centers = cluster_info$best_k, nstart = 10)
  cluster_labels <- final_model$cluster
  }
  feature_desc <- if (grepl("observed", output_prefix)) "#sbst/#ori" else "log2((#sbst/#ori)/mean)"
  if (length(plot_levels) > 0) {
    for (level in plot_levels) {
      color_labels <- vapply(metadata_records[indices], function(rec) rec$classifications[[level]] %||% "Unknown", character(1))
      plot_path <- file.path(output_dir, sprintf("pca_2d_clusters_%s_%s.png", output_prefix, gsub(" ", "_", level)))
      plot_pca_scatter(plot_path, pca_data, species_labels[indices], color_labels, cluster_labels, level, feature_desc)
    }
  }
  cluster_color_labels <- as.character(cluster_labels)
  cluster_plot_path <- file.path(output_dir, sprintf("pca_2d_clusters_%s_cluster.png", output_prefix))
  plot_pca_scatter(cluster_plot_path, pca_data, species_labels[indices], cluster_color_labels, cluster_labels, "cluster", feature_desc)
  save_results(
    file.path(output_dir, sprintf("clustering_results_%s.csv", output_prefix)),
    files[indices],
    pca_data,
    cluster_labels,
    metadata_records[indices],
    classification_levels
  )
  list(best_k = cluster_info$best_k, processed = length(indices))
}

process_variant <- function(vectors, feature_order, indices_all, indices_filtered, variant_label, plot_levels, metadata_records, files, classification_levels, species_labels, max_clusters, random_state, output_dir) {
  normal_prefix <- variant_label
  filtered_prefix <- sprintf("%s_filtered", variant_label)
  normal_result <- run_pca_for_indices(
    vectors, feature_order, indices_all, normal_prefix, plot_levels, metadata_records, files, classification_levels, species_labels, max_clusters, random_state, output_dir
  )
  filtered_result <- list(best_k = NA, processed = 0, message = NULL)
  if (length(indices_filtered) >= 2) {
    filtered_inner <- run_pca_for_indices(
      vectors, feature_order, indices_filtered, filtered_prefix, plot_levels, metadata_records, files, classification_levels, species_labels, max_clusters, random_state, output_dir
    )
    filtered_result$best_k <- filtered_inner$best_k
    filtered_result$processed <- filtered_inner$processed
  } else {
    message_text <- if (length(indices_filtered) == 0) {
      sprintf("%s filtered PCA skipped: no samples remained after applying the filter.", variant_label)
    } else {
      sprintf("%s filtered PCA requires at least two samples; skipping filtered computation.", variant_label)
    }
    filtered_result$message <- message_text
    write_placeholder_plot(file.path(output_dir, sprintf("clustering_metrics_%s.png", filtered_prefix)), message_text)
    write_empty_results_csv(file.path(output_dir, sprintf("clustering_results_%s.csv", filtered_prefix)), classification_levels)
    if (length(plot_levels) > 0) {
      for (level in plot_levels) {
        placeholder_path <- file.path(output_dir, sprintf("pca_2d_clusters_%s_%s.png", filtered_prefix, gsub(" ", "_", level)))
        write_placeholder_plot(placeholder_path, message_text)
      }
    }
  }
  list(
    variant = variant_label,
    normal = normal_result,
    filtered = filtered_result
  )
}

main <- function() {
  args <- parse_args()
  input_root <- normalizePath(args$input_dir, mustWork = TRUE)
  if (!dir.exists(input_root)) stop(sprintf("Input directory does not exist: %s", input_root))
  if (!is.null(args$output_dir)) {
    output_dir <- normalizePath(args$output_dir, mustWork = FALSE)
  } else {
    output_dir <- normalizePath(file.path(input_root, "pca"), mustWork = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

  dataset_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
  dataset_dirs <- dataset_dirs[!grepl("^[0-9]+$", basename(dataset_dirs))]
  if (length(dataset_dirs) == 0) dataset_dirs <- input_root

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
    manifest_bundle <- tryCatch(load_manifest_lookup(manifest_path), error = function(e) {
      warn(sprintf("%s: Failed to load manifest %s: %s", dataset_name, manifest_path, e$message))
      NULL
    })
    if (is.null(manifest_bundle)) next
    slot_map <- manifest_bundle$slot_map
    filter_info <- evaluate_dataset_filter(slot_map, latest_run)
    for (issue in filter_info$issues) warn(sprintf("%s: %s", dataset_name, issue))
    dataset_status[[dataset_name]] <- filter_info$filter_status
    dataset_filter_flag <- isTRUE(filter_info$filter_status$excluded)

    tsv_files <- Sys.glob(file.path(latest_run, args$tsv_pattern))
    if (length(tsv_files) == 0) {
      warn(sprintf("%s: No TSV files matching '%s' in %s.", dataset_name, args$tsv_pattern, latest_run))
      next
    }
    samples_before <- length(logratio_vectors)
    for (tsv_file in sort(tsv_files)) {
      parts <- strsplit(basename(tsv_file), "_")[[1]]
      if (length(parts) < 4) {
        warn(sprintf("%s: Unexpected TSV filename format: %s", dataset_name, basename(tsv_file)))
        next
      }
      accession <- paste(parts[1:2], collapse = "_")
      short_name <- parts[3]
      date <- parts[4]
      key <- paste(accession, short_name, sep = "::")
      manifest_entry <- manifest_bundle$lookup[[key]]
      if (is.null(manifest_entry)) {
        warn(sprintf("%s: No manifest entry for accession %s / short_name %s (file %s).", dataset_name, accession, short_name, basename(tsv_file)))
        next
      }
      metadata_path <- manifest_entry$metadata_path
      metadata_json <- tryCatch(load_json_cached(metadata_path, metadata_cache), error = function(e) {
        warn(sprintf("%s: Failed to load metadata JSON %s: %s", dataset_name, metadata_path, e$message))
        NULL
      })
      if (is.null(metadata_json)) next
      taxonomy_path <- manifest_entry$taxonomy_path
      ensure_taxonomy_json(taxonomy_path, metadata_json, short_name, accession)
      taxonomy_json <- tryCatch(load_json_cached(taxonomy_path, taxonomy_cache), error = function(e) {
        warn(sprintf("%s: Failed to load taxonomy JSON %s: %s", dataset_name, taxonomy_path, e$message))
        list()
      })
      genus_info <- extract_genus(metadata_json, short_name, accession)
      classification_map <- extract_classification_map(taxonomy_json, sprintf("%s (%s)", short_name, accession))
      observed_levels <- unique(c(observed_levels, names(classification_map)))

      vectors <- tryCatch(load_tsv_vectors(tsv_file, reference_order), error = function(e) {
        warn(sprintf("%s: Failed to parse TSV %s: %s", dataset_name, basename(tsv_file), e$message))
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
        dataset = dataset_name,
        filter_excluded = dataset_filter_flag
    )
  }
    added <- length(logratio_vectors) - samples_before
    if (added == 0) warn(sprintf("%s: No usable TSV files after validation; dataset skipped.", dataset_name))
  }

  if (length(logratio_vectors) < 2) stop("At least two valid TSV samples are required for PCA.")

  classification_levels <- intersect(CLASSIFICATION_ORDER, unique(observed_levels))
  remaining_levels <- setdiff(unique(observed_levels), classification_levels)
  classification_levels <- c(classification_levels, sort(remaining_levels))
  if (length(classification_levels) == 0) classification_levels <- "taxonomy"
  plot_levels <- setdiff(classification_levels, "species")

  indices_all <- seq_along(logratio_vectors)
  indices_filtered <- which(!vapply(metadata_records, function(rec) rec$filter_excluded, logical(1)))

  logratio_result <- process_variant(
    logratio_vectors,
    reference_order,
    indices_all,
    indices_filtered,
    "logratio",
    plot_levels,
    metadata_records,
    files,
    classification_levels,
    species_labels,
    args$max_clusters,
    args$random_state,
    output_dir
  )

  observed_result <- process_variant(
    observed_vectors,
    reference_order,
    indices_all,
    indices_filtered,
    "observed",
    plot_levels,
    metadata_records,
    files,
    classification_levels,
    species_labels,
    args$max_clusters,
    args$random_state,
    output_dir
  )

  message("PCA clustering completed successfully.")
  message(sprintf("Datasets processed: %d", length(dataset_status)))
  message(sprintf("Total samples processed: %d", length(indices_all)))
  for (result in list(logratio_result, observed_result)) {
    message(sprintf(
      "%s variant — processed: %d, optimal k: %s",
      stringr::str_to_title(result$variant),
      result$normal$processed,
      result$normal$best_k
    ))
    if (result$filtered$processed > 0) {
      message(sprintf("  Filtered processed files: %d, optimal k: %s", result$filtered$processed, result$filtered$best_k))
    } else if (!is.null(result$filtered$message)) {
      message(sprintf("  %s", result$filtered$message))
    }
  }
  if (length(plot_levels) == 0) {
    message("No classification plots generated (species-only classifications).")
  } else {
    message(sprintf("Generated classification plots for levels: %s", paste(plot_levels, collapse = ", ")))
  }
  if (length(dataset_status) > 0) {
    message("Filter evaluation summary:")
    for (dataset_name in names(dataset_status)) {
      status <- dataset_status[[dataset_name]]
      message(sprintf(
        "  %s: excluded=%s genus_condition=%s identity_condition=%s",
        dataset_name,
        status$excluded,
        status$genus_condition_met,
        status$identity_condition_met
      ))
    }
  }
  message(sprintf("Results saved to: %s", output_dir))
}

tryCatch(
  main(),
  error = function(e) {
    message(sprintf("Error: %s", e$message))
    quit(status = 1)
  }
)
