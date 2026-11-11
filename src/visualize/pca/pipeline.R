# Coordinate computation, plotting, and failure-safe run reporting.

run_pca_for_indices <- function(vectors, feature_order, indices, output_prefix, plot_levels, metadata_records, files, classification_levels, species_labels, max_clusters, random_state, output_dir, cluster_boundary) {
  if (length(indices) < 2) stop(sprintf("%s PCA requires at least two samples.", output_prefix))
  projection <- compute_pca_projection(vectors, feature_order, indices)
  pca <- projection$pca
  pca_data <- projection$data
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

      # Original: Labels + Boundary (if enabled)
      plot_path <- file.path(output_dir, sprintf("pca_2d_clusters_%s_%s.png", output_prefix, gsub(" ", "_", level)))
      plot_pca_scatter(plot_path, pca_data, species_labels[indices], color_labels, cluster_labels, level, feature_desc, cluster_boundary = cluster_boundary, show_labels = TRUE)

      # Variant 1: No Labels + Boundary (if enabled) (_noName)
      plot_path_noName <- file.path(output_dir, sprintf("pca_2d_clusters_%s_%s_noName.png", output_prefix, gsub(" ", "_", level)))
      plot_pca_scatter(plot_path_noName, pca_data, species_labels[indices], color_labels, cluster_labels, level, feature_desc, cluster_boundary = cluster_boundary, show_labels = FALSE)

      # Variant 2: No Labels + No Boundary (_pcaOnly_noName)
      plot_path_pcaOnly <- file.path(output_dir, sprintf("pca_2d_clusters_%s_%s_pcaOnly_noName.png", output_prefix, gsub(" ", "_", level)))
      plot_pca_scatter(plot_path_pcaOnly, pca_data, species_labels[indices], color_labels, cluster_labels, level, feature_desc, cluster_boundary = "none", show_labels = FALSE)
    }
  }
  cluster_color_labels <- as.character(cluster_labels)
  cluster_plot_path <- file.path(output_dir, sprintf("pca_2d_clusters_%s_cluster.png", output_prefix))
  plot_pca_scatter(cluster_plot_path, pca_data, species_labels[indices], cluster_color_labels, cluster_labels, "cluster", feature_desc, cluster_boundary = "none")
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

process_variant <- function(vectors, feature_order, indices_all, variant_label, plot_levels, metadata_records, files, classification_levels, species_labels, max_clusters, random_state, output_dir, cluster_boundary) {
  result <- run_pca_for_indices(vectors, feature_order, indices_all, variant_label, plot_levels, metadata_records, files, classification_levels, species_labels, max_clusters, random_state, output_dir, cluster_boundary)
  list(variant=variant_label, normal=result)
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

  if (length(list.files(output_dir, pattern = "_filtered"))) {
    warn("Existing _filtered files are artifacts from an earlier run and are not updated.")
  }
  option_keys <- c(
    "idt_threshold", "ingroup_idt_threshold", "outgroup_idt_threshold",
    "min_idt_gap_between_out_in"
  )
  options <- args[option_keys]
  effective <- list(
    ingroup_idt_threshold = options$ingroup_idt_threshold %||% options$idt_threshold,
    outgroup_idt_threshold = options$outgroup_idt_threshold %||% options$idt_threshold,
    min_idt_gap_between_out_in = options$min_idt_gap_between_out_in
  )
  settings <- list(
    arguments = commandArgs(trailingOnly = TRUE),
    requested_thresholds = options,
    effective_thresholds = effective,
    disabled_conditions = names(effective)[vapply(effective, is.null, logical(1))],
    identity_source = tail(system2("python3", c(
      shQuote(file.path(dirname(.pca_manifest_helper), "train_identity.py")), "--source"
    ), stdout = TRUE), 1),
    comparison = ">",
    identity_unit = "percent",
    figure_background = "transparent",
    gap_unit = "percentage points",
    genus_waiver = "same-genus ingroups and different-genus outgroup waive only gap",
    input_dir = input_root,
    output_dir = output_dir,
    pca_options = args,
    started_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    status = "running"
  )
  audit <- new_pca_audit(settings, output_dir)
  on.exit(save_pca_audit(audit), add = TRUE)
  save_pca_audit(audit)
  audit$settings$status <- "failed"

  inputs <- collect_pca_inputs(input_root, args, options, audit)

  if (length(inputs$logratio_vectors) < 2) stop("At least two valid TSV samples are required for PCA.")

  classification_levels <- intersect(CLASSIFICATION_ORDER, unique(inputs$observed_levels))
  remaining_levels <- setdiff(unique(inputs$observed_levels), classification_levels)
  classification_levels <- c(classification_levels, sort(remaining_levels))
  if (length(classification_levels) == 0) classification_levels <- "taxonomy"
  plot_levels <- setdiff(classification_levels, "species")

  indices_all <- seq_along(inputs$logratio_vectors)

  logratio_result <- process_variant(
    inputs$logratio_vectors,
    inputs$reference_order,
    indices_all,
    "logratio",
    plot_levels,
    inputs$metadata_records,
    inputs$files,
    classification_levels,
    inputs$species_labels,
    args$max_clusters,
    args$random_state,
    output_dir,
    args$cluster_boundary
  )

  observed_result <- process_variant(
    inputs$observed_vectors,
    inputs$reference_order,
    indices_all,
    "observed",
    plot_levels,
    inputs$metadata_records,
    inputs$files,
    classification_levels,
    inputs$species_labels,
    args$max_clusters,
    args$random_state,
    output_dir,
    args$cluster_boundary
  )

  message("PCA clustering completed successfully.")
  message(sprintf("Datasets processed: %d", length(inputs$dataset_status)))
  message(sprintf("Total samples processed: %d", length(indices_all)))
  for (result in list(logratio_result, observed_result)) {
    message(sprintf(
      "%s variant — processed: %d, optimal k: %s",
      stringr::str_to_title(result$variant),
      result$normal$processed,
      result$normal$best_k
    ))

  }
  if (length(plot_levels) == 0) {
    message("No classification plots generated (species-only classifications).")
  } else {
    message(sprintf("Generated classification plots for levels: %s", paste(plot_levels, collapse = ", ")))
  }
  audit$settings$status <- "completed"
  message(sprintf("Results saved to: %s", output_dir))
}
