#!/usr/bin/env Rscript
.shape_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  if (length(files)) tail(files, 1)[[1]] else sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
.shape_dir <- dirname(normalizePath(.shape_file))
source(file.path(.shape_dir, "transparent_background.R"), local = TRUE)
for (part in c("core", "io", "legacy", "mixed", "plots")) source(file.path(.shape_dir, "spectrum_shape", paste0(part, ".R")))
shape_main <- function() {
  parser <- argparse::ArgumentParser(description = "Cluster 96-channel log2-ratio spectra and annotate taxonomy.")
  parser$add_argument("input_dir", help = "Root containing v2 runs or legacy trio directories")
  parser$add_argument("--output-dir", default = NULL, help = "Fresh empty output directory; default INPUT_DIR/spectrum_shape")
  parser$add_argument("--input-mode", default = "v2", choices = c("v2", "legacy", "mixed"))
  parser$add_argument("--run-pattern", default = "*")
  parser$add_argument("--region", default = NULL, choices = c("whole", "noncoding"))
  parser$add_argument("--tsv-pattern", default = NULL, choices = c("*_ncds.tsv", "*.tsv"))
  parser$add_argument("--taxonomy-levels", default = "phylum,class,order", help = "Comma-separated taxonomy ranks for annotation only")
  parser$add_argument("--clusters", default = "auto", help = "auto: maximize silhouette, or a positive integer cluster count")
  parser$add_argument("--max-clusters", default = 10L, type = "integer")
  parser$add_argument("--random-state", default = 42L, type = "integer")
  parser$add_argument("--method", default = "ward", choices = c("ward", "correlation"))
  args <- parser$parse_args()
  args$region <- shape_region(args$region, args$tsv_pattern)
  args$tsv_pattern <- if (args$region == "whole") "*.tsv" else "*_ncds.tsv"
  if (args$clusters != "auto" && (!grepl("^[1-9][0-9]*$", args$clusters) || is.na(suppressWarnings(as.integer(args$clusters))))) stop("Invalid --clusters")
  if (args$max_clusters < 2) stop("--max-clusters must be at least 2")
  levels <- unique(trimws(strsplit(args$taxonomy_levels, ",", fixed = TRUE)[[1]]))
  if (!length(levels) || any(!levels %in% c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"))) stop("Invalid taxonomy levels")
  # Prefix taxonomy columns to keep species metadata distinct from species taxonomy.
  input <- normalizePath(args$input_dir, mustWork = TRUE)
  out <- if (is.null(args$output_dir)) file.path(input, "spectrum_shape") else args$output_dir
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(out)
  if (length(existing)) stop("Output directory must be empty to prevent stale results: ", out)
  data <- if (args$input_mode == "mixed") {
    shape_mixed_ingest(input, args$region, levels, args$method, args$run_pattern, file.path(.shape_dir, "..", "lib", "run_manifest.py"))
  } else if (args$input_mode == "legacy") {
    shape_legacy_ingest(input, args$region, levels, args$method, args$run_pattern)
  } else {
    shape_ingest(input, args$tsv_pattern, levels, file.path(.shape_dir, "..", "lib", "run_manifest.py"), args$method)
  }
  if (!is.null(data$run_decisions)) shape_csv(data$run_decisions, file.path(out, "run_decisions.csv"))
  shape_csv(data$records, file.path(out, "samples.csv"))
  shape_csv(data$records[!data$records$included, , drop = FALSE], file.path(out, "exclusions.csv"))
  counts <- do.call(rbind, lapply(levels, function(level) {
    if (!nrow(data$records)) {
      return(data.frame(level = character(), taxon = character(), included = logical(), n = integer()))
    }
    z <- aggregate(rep(1, nrow(data$records)), list(taxon = data$records[[if (level == "species") "taxonomy_species" else level]], included = data$records$included), sum)
    names(z)[3] <- "n"
    z$level <- level
    z
  }))
  shape_csv(counts, file.path(out, "taxonomy_counts.csv"))
  x <- data$matrix
  shape_matrix_csv(x, file.path(out, "log2_matrix.csv"))
  jsonlite::write_json(list(arguments = args, input = input, runs = data$runs, provenance = switch(args$input_mode, v2 = "v2 manifest validation", legacy = "legacy input; v2 completion and provenance not verified", mixed = "mixed input; selected v2 runs validated, legacy runs not v2-verified"), formula = "log2((mutNum/totalRootNum)/mean(positive rates)); no pseudocount", seed = args$random_state, figure_background = "transparent"), file.path(out, "settings.json"), auto_unbox = TRUE, pretty = TRUE)
  if (nrow(x) < 2) stop("At least two valid samples are required; see exclusions.csv")
  packages <- c("ComplexHeatmap", "circlize", "ggplot2", "uwot", "cluster", "argparse", "jsonlite")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  shape_csv(data.frame(package = packages, version = vapply(packages, function(p) as.character(packageVersion(p)), character(1))), file.path(out, "versions.csv"))
  result <- shape_cluster(x, args$method, args$clusters, args$max_clusters)
  shape_matrix_csv(as.matrix(result$distance), file.path(out, "distances.csv"))
  shape_csv(data.frame(sample_id = rownames(x), cluster = result$membership), file.path(out, "clusters.csv"))
  shape_csv(result$scores, file.path(out, "cluster_scores.csv"))
  shape_csv(data.frame(cluster = seq_along(result$representatives), sample_id = result$representatives), file.path(out, "representatives.csv"))
  meta <- data$records[match(rownames(x), data$records$sample_id), , drop = FALSE]
  composition <- do.call(rbind, lapply(levels, function(level) {
    z <- aggregate(rep(1, nrow(meta)), list(cluster = result$membership, taxon = meta[[if (level == "species") "taxonomy_species" else level]]), sum)
    names(z)[3] <- "n"
    z$level <- level
    z$fraction <- z$n / as.numeric(table(result$membership)[as.character(z$cluster)])
    z
  }))
  shape_csv(composition, file.path(out, "cluster_taxonomy.csv"))
  embedding <- shape_umap(x, args$method, args$random_state)
  if (is.null(embedding)) {
    writeLines("UMAP skipped: fewer than four samples", file.path(out, "umap_status.txt"))
  } else {
    shape_csv(data.frame(sample_id = rownames(x), UMAP1 = embedding[, 1], UMAP2 = embedding[, 2]), file.path(out, "umap.csv"))
  }
  writeLines(if (result$k == 1 && args$clusters == "auto") "Auto silhouette unavailable: two samples or all distances zero; one cluster" else paste("Selected clusters:", result$k), file.path(out, "cluster_status.txt"))
  shape_plots(x, meta, result, embedding, levels, out)
  message("Spectrum shape results: ", normalizePath(out))
}
if (sys.nframe() == 0) {
  tryCatch(shape_main(), error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1)
  })
}
