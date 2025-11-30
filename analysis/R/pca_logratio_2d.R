#!/usr/bin/env Rscript

# PCA clustering pipeline that mirrors scripts/analysis/pca_logratio_2d.py.
# Usage example:
#   Rscript scripts/analysis/R/pca_logratio_2d.R \
#     --tsv-dir /path/to/tsv_dir --output-dir /path/to/output_dir
# Required packages: argparse, jsonlite, tidyverse, cluster, scales, factoextra, ggpubr.
# The script calls `datasets summary genome accession "<accession>"` when metadata JSON is missing.

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
  library(factoextra)
  library(ggpubr)
  library(rlang)
})

EPS <- 1e-12

# Suppress R CMD check notes for dplyr/ggplot2 pronouns
utils::globalVariables(c("mutNum", "totalRootNum", "mutType", "k", "PC1", "PC2", "genus", "species"))

parse_args <- function() {
  parser <- ArgumentParser(description = "Run 2D PCA clustering on substitution TSV files (logRatio vectors).")
  parser$add_argument("--tsv-dir", required = TRUE, help = "Directory containing substitution TSV files.")
  parser$add_argument("--output-dir", required = FALSE, help = "Directory to write PCA outputs (default: <tsv-dir>/pca_logratio_output).")
  parser$add_argument("--tsv-pattern", default = "*.tsv", help = "Glob pattern for TSV files (default: *.tsv).")
  parser$add_argument("--max-clusters", type = "integer", default = 10L, help = "Maximum number of clusters to evaluate.")
  parser$add_argument("--random-state", type = "integer", default = 42L, help = "Random seed for clustering.")
  args <- parser$parse_args()
  args
}

collect_tsv_files <- function(tsv_dir, pattern) {
  files <- list.files(path = tsv_dir, pattern = glob2rx(pattern), full.names = TRUE, recursive = FALSE)
  if (length(files) == 0) {
    stop(sprintf("No TSV files matching pattern '%s' were found in %s", pattern, tsv_dir))
  }
  sort(files)
}

parse_filename_parts <- function(file_path) {
  name_parts <- strsplit(sub("\\.tsv$", "", basename(file_path)), "_")[[1]]
  if (length(name_parts) < 4) {
    stop(sprintf("Unexpected filename format for TSV file: %s", basename(file_path)))
  }
  list(
    accession = paste(name_parts[1:2], collapse = "_"),
    short_name = name_parts[3],
    date = name_parts[4]
  )
}

run_datasets_cli <- function(accession) {
  command <- sprintf('datasets summary genome accession "%s"', accession)
  result <- system(command, intern = TRUE, ignore.stderr = TRUE)
  if (length(result) == 0) {
    stop(sprintf("Empty metadata response for accession %s", accession))
  }
  fromJSON(paste(result, collapse = "\n"))
}

ensure_metadata <- function(tsv_path, accession) {
  metadata_path <- file.path(dirname(tsv_path), paste0(accession, ".json"))
  if (file.exists(metadata_path)) {
    return(fromJSON(metadata_path, simplifyVector = FALSE))
  }
  metadata <- run_datasets_cli(accession)
  write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE)
  metadata
}

extract_genus_label <- function(metadata, short_name, accession) {
  reports <- metadata$reports
  if (!is.null(reports) && length(reports) > 0) {
    organism <- reports[[1]]$organism
    if (!is.null(organism$organism_name)) {
      species_name <- organism$organism_name
      genus <- strsplit(species_name, " +")[[1]][1]
      return(list(genus = genus, species = species_name))
    }
  }
  list(genus = "Unknown", species = accession)
}

compute_feature_vectors <- function(df) {
  if (!all(c("mutType", "mutNum", "totalRootNum") %in% names(df))) {
    stop("TSV file is missing required columns (mutType, mutNum, totalRootNum)")
  }
  df <- df %>%
    mutate(
      mutNum = as.numeric(.data$mutNum),
      totalRootNum = as.numeric(.data$totalRootNum),
      observed = if_else(.data$totalRootNum != 0, .data$mutNum / .data$totalRootNum, 0)
    )
  mean_observed <- ifelse(nrow(df) > 0, mean(df$observed, na.rm = TRUE), 0)
  adjusted_mean <- max(mean_observed, EPS)
  ratios <- (df$observed + EPS) / adjusted_mean
  list(
    logratio = log2(ratios),
    observed = df$observed
  )
}

load_tsv_vectors <- function(file_path, reference_order = NULL) {
  df <- read_tsv(file_path, show_col_types = FALSE)
  if (is.null(reference_order)) {
    df_sorted <- df %>% arrange(.data$mutType)
    vectors <- compute_feature_vectors(df_sorted)
    order <- df_sorted$mutType
    return(list(logratio = vectors$logratio, observed = vectors$observed, order = order))
  }
  df_indexed <- df %>%
    filter(.data$mutType %in% reference_order) %>%
    arrange(match(.data$mutType, reference_order))
  if (!all(reference_order %in% df_indexed$mutType)) {
    stop(sprintf("TSV file %s does not contain expected mutType values", basename(file_path)))
  }
  vectors <- compute_feature_vectors(df_indexed)
  list(logratio = vectors$logratio, observed = vectors$observed, order = reference_order)
}

determine_clusters <- function(data, max_clusters, seed) {
  sample_count <- nrow(data)
  if (sample_count < 2) stop("At least two samples are required for clustering.")
  adjusted_max <- min(max_clusters, sample_count)
  if (adjusted_max < 2) adjusted_max <- 2

  set.seed(seed)
  inertia_values <- numeric()
  silhouette_scores <- numeric()

  for (k in 2:adjusted_max) {
    model <- kmeans(data, centers = k, nstart = 10)
    inertia_values <- c(inertia_values, model$tot.withinss)
    sil <- silhouette(model$cluster, dist(data))
    silhouette_scores <- c(silhouette_scores, mean(sil[, "sil_width"]))
  }

  best_k <- which.max(silhouette_scores) + 1
  list(best_k = best_k, inertia = inertia_values, silhouette = silhouette_scores)
}

plot_metrics <- function(output_dir, inertia, silhouette, suffix) {
  ks <- seq_along(inertia) + 1
  metrics_df <- tibble(
    k = ks,
    inertia = inertia,
    silhouette = silhouette
  )

  p1 <- ggplot(metrics_df, aes(x = .data$k, y = .data$inertia)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    labs(title = "Elbow Method", x = "Number of Clusters", y = "Inertia") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "#d0d0d0", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )

  p2 <- ggplot(metrics_df, aes(x = .data$k, y = .data$silhouette)) +
    geom_line(color = "red") +
    geom_point(color = "red") +
    labs(title = "Silhouette Scores", x = "Number of Clusters", y = "Silhouette Score") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "#d0d0d0", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )

  file_name <- if (suffix == "logratio") "clustering_metrics.png" else sprintf("clustering_metrics_%s.png", suffix)

  ggsave(
    file.path(output_dir, file_name),
    ggpubr::ggarrange(p1, p2, nrow = 1),
    dpi = 300,
    width = 12,
    height = 5
  )
}

plot_pca_scatter <- function(output_dir, pca_data, labels, genus_labels, cluster_labels, suffix) {
  df <- tibble(
    PC1 = pca_data[, 1],
    PC2 = pca_data[, 2],
    species = stringr::str_replace_all(labels, "_", " "),
    genus = genus_labels,
    cluster = factor(cluster_labels)
  )

  palette <- hue_pal()(length(unique(df$genus)))
  names(palette) <- unique(df$genus)

  subtitle_text <- if (suffix == "logratio") "Features: log2((#sbst/#ori)/mean of #sbst/#ori)" else "Features: #sbst/#ori"

  p <- ggplot(df, aes(x = .data$PC1, y = .data$PC2, color = .data$genus)) +
    geom_point(size = 3, alpha = 0.8) +
    ggrepel::geom_text_repel(
      aes(label = .data$species),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = Inf,
      segment.color = "gray50",
      segment.alpha = 0.7,
      fontface = "italic",
      force = 3,
      min.segment.length = 0
    ) +
    scale_color_manual(values = palette) +
    labs(
      title = sprintf("K-means Clustering (k=%d)", length(unique(cluster_labels))),
      subtitle = subtitle_text,
      x = "PC1",
      y = "PC2"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "#d0d0d0", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )

  file_name <- if (suffix == "logratio") "pca_2d_clusters.png" else sprintf("pca_2d_clusters_%s.png", suffix)

  ggsave(file.path(output_dir, file_name), p, dpi = 300, width = 10, height = 8)
}

run_pca_workflow <- function(data_matrix, feature_suffix, output_dir, species_labels, genus_labels, max_clusters, random_state) {
  scaled_data <- scale(data_matrix)
  pca <- prcomp(scaled_data, center = FALSE, scale. = FALSE)
  pca_data <- pca$x[, 1:2, drop = FALSE]

  cluster_info <- determine_clusters(pca_data, max_clusters, random_state)
  plot_metrics(output_dir, cluster_info$inertia, cluster_info$silhouette, feature_suffix)

  set.seed(random_state)
  final_model <- kmeans(pca_data, centers = cluster_info$best_k, nstart = 10)
  cluster_labels <- final_model$cluster

  plot_pca_scatter(output_dir, pca_data, species_labels, genus_labels, cluster_labels, feature_suffix)

  list(
    cluster_labels = cluster_labels,
    pca_data = pca_data,
    best_k = cluster_info$best_k
  )
}

save_results <- function(output_dir, files, metadata_records, logratio_results, observed_results) {
  df <- tibble(
    file = basename(files),
    accession = map_chr(metadata_records, "accession"),
    short_name = map_chr(metadata_records, "short_name"),
    date = map_chr(metadata_records, "date"),
    genus = map_chr(metadata_records, "genus"),
    species = map_chr(metadata_records, "species"),
    cluster_logratio = as.integer(logratio_results$cluster_labels),
    cluster_observed = as.integer(observed_results$cluster_labels),
    PC1_logratio = logratio_results$pca_data[, 1],
    PC2_logratio = logratio_results$pca_data[, 2],
    PC1_observed = observed_results$pca_data[, 1],
    PC2_observed = observed_results$pca_data[, 2]
  )
  write_csv(df, file.path(output_dir, "clustering_results.csv"))
}

run <- function() {
  args <- parse_args()

  tsv_dir <- normalizePath(args$tsv_dir)
  if (!dir.exists(tsv_dir)) {
    stop(sprintf("TSV directory does not exist: %s", tsv_dir))
  }
  output_dir <- if (!is.null(args$output_dir)) normalizePath(args$output_dir) else file.path(tsv_dir, "pca_logratio_output")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  files <- collect_tsv_files(tsv_dir, args$tsv_pattern)

  logratio_vectors <- list()
  observed_vectors <- list()
  labels <- character()
  genus_labels <- character()
  metadata_records <- list()
  reference_order <- NULL

  for (file in files) {
    parts <- parse_filename_parts(file)
    metadata <- ensure_metadata(file, parts$accession)
    label_info <- extract_genus_label(metadata, parts$short_name, parts$accession)

    loaded <- load_tsv_vectors(file, reference_order)
    reference_order <- loaded$order

    logratio_vectors[[length(logratio_vectors) + 1]] <- loaded$logratio
    observed_vectors[[length(observed_vectors) + 1]] <- loaded$observed
    labels <- c(labels, label_info$species)
    genus_labels <- c(genus_labels, label_info$genus)
    metadata_records[[length(metadata_records) + 1]] <- list(
      accession = parts$accession,
      short_name = parts$short_name,
      date = parts$date,
      genus = label_info$genus,
      species = label_info$species
    )
  }

  logratio_matrix <- do.call(rbind, logratio_vectors)
  observed_matrix <- do.call(rbind, observed_vectors)

  logratio_results <- run_pca_workflow(
    logratio_matrix,
    "logratio",
    output_dir,
    labels,
    genus_labels,
    args$max_clusters,
    args$random_state
  )

  observed_results <- run_pca_workflow(
    observed_matrix,
    "observed",
    output_dir,
    labels,
    genus_labels,
    args$max_clusters,
    args$random_state
  )

  save_results(output_dir, files, metadata_records, logratio_results, observed_results)

  message("PCA clustering completed successfully.")
  message(sprintf("Processed files: %d", length(files)))
  message(sprintf("Optimal clusters (logratio k): %d", logratio_results$best_k))
  message(sprintf("Optimal clusters (observed k): %d", observed_results$best_k))
  message(sprintf("Results saved to: %s", output_dir))
}

tryCatch(
  run(),
  error = function(e) {
    message(sprintf("Error: %s", e$message))
    quit(status = 1)
  }
)
