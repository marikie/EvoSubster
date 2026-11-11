#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))

test_modules <- function() {
  work <- tempfile("pca-modules-", tmpdir = file.path(root, "tmp", "codex"))
  dir.create(work)
  previous_dir <- getwd()
  on.exit(setwd(previous_dir), add = TRUE)
  setwd(work)

  # Source from another directory without running main or leaking its helpers.
  pca <- new.env(parent = globalenv())
  source(file.path(root, "src", "visualize", "pca_2d.R"), local = pca)
  stopifnot(
    identical(environment(pca$main), pca),
    identical(environment(pca$collect_pca_inputs), pca),
    identical(environment(pca$save_pca_audit), pca),
    identical(pca$.pca_manifest_helper, file.path(root, "src", "lib", "pca_input_adapter.py")),
    length(list.files(work)) == 0
  )

  # Capture plots to check geometry without raster or label-placement variation.
  pca$ggsave <- function(filename, plot, ...) plot
  coordinates <- matrix(c(0, 0, 2, 0, 0, 2, 10, 10, 12, 10, 10, 12), ncol = 2, byrow = TRUE)
  labels <- rep(c("Mucoromycota", "Apicomplexa"), each = 3)
  clusters <- rep(1:2, each = 3)
  build_plot <- function(boundary) {
    pca$plot_pca_scatter("unused.png", coordinates, letters[1:6], labels, clusters,
                         "phylum", "test", cluster_boundary = boundary, show_labels = FALSE)
  }
  plain <- ggplot2::ggplot_build(build_plot("none"))
  stopifnot(length(plain$data) == 1, identical(plain$data[[1]]$colour, rep(c("#56B4E9", "#E69F00"), each = 3)))
  hull <- ggplot2::ggplot_build(build_plot("hull"))
  stopifnot(length(hull$data) == 2)
  hull_data <- hull$data[[2]]
  for (group in 1:2) {
    actual <- unique(hull_data[hull_data$group == group, c("x", "y")])
    offset <- (group - 1) * 10
    expected <- data.frame(x = c(-1 / 15, 32 / 15, -1 / 15) + offset,
                           y = c(-1 / 15, -1 / 15, 32 / 15) + offset)
    actual <- actual[order(actual$x, actual$y), ]
    expected <- expected[order(expected$x, expected$y), ]
    stopifnot(isTRUE(all.equal(unname(as.matrix(actual)), unname(as.matrix(expected)), tolerance = 1e-12)))
  }
  ellipse_coordinates <- rbind(coordinates[1:3, ], c(2, 2), coordinates[4:6, ], c(12, 12))
  ellipse_plot <- pca$plot_pca_scatter(
    "unused.png", ellipse_coordinates, letters[1:8], rep(c("Mucoromycota", "Apicomplexa"), each = 4),
    rep(1:2, each = 4), "phylum", "test", cluster_boundary = "ellipse", show_labels = FALSE
  )
  ellipse <- ggplot2::ggplot_build(ellipse_plot)
  stopifnot(length(ellipse$data) == 2, length(unique(ellipse$data[[2]]$group)) == 2,
            all(is.finite(ellipse$data[[2]]$x)), all(is.finite(ellipse$data[[2]]$y)))

  setwd(previous_dir)
  unlink(work, recursive = TRUE)
}

test_modules()
cat("All PCA module and boundary tests passed.\n")
