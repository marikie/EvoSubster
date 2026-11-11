# Draw the standard PCA scatter and cluster-metric figures.

plot_metrics <- function(metrics_path, inertia, silhouette) {
  if (length(inertia) == 0 && length(silhouette) == 0) {
    ggsave(metrics_path, ggplot() +
      theme_void() + transparent_background_theme(), width = 8, height = 4, dpi = 300, bg = "transparent")
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
  p1 <- p1 + transparent_background_theme()
  p2 <- p2 + transparent_background_theme()
  # Build before ggsave opens its device, preserving the original text metrics.
  arranged <- gridExtra::arrangeGrob(p1, p2, nrow = 1)
  ggsave(metrics_path, arranged, dpi = 300, width = 12, height = 5, bg = "transparent")
}

build_classification_palette <- function(labels, label_name) {
  labels <- unique(as.character(labels))
  palette <- scales::hue_pal()(length(labels))
  names(palette) <- labels
  if (tolower(label_name %||% "") != "phylum") {
    return(palette)
  }

  known_labels <- intersect(labels, names(PHYLUM_COLORS))
  palette[known_labels] <- PHYLUM_COLORS[known_labels]
  palette
}

plot_pca_scatter <- function(output_path, pca_data, species_labels, color_labels, cluster_labels, label_name, feature_desc, cluster_boundary = "none", show_labels = TRUE) {
  df <- tibble(
    PC1 = pca_data[, 1],
    PC2 = pca_data[, 2],
    species = stringr::str_replace_all(species_labels, "_", " "),
    label = color_labels,
    cluster = factor(cluster_labels)
  )
  palette <- build_classification_palette(df$label, label_name)
  classification_label <- stringr::str_to_title(label_name)
  boundary_mode <- tolower(cluster_boundary %||% "none")
  p <- ggplot(df, aes(x = .data$PC1, y = .data$PC2, color = .data$label)) +
    geom_point(size = 3, alpha = 0.85) +
    {
      if (boundary_mode == "ellipse") {
        df_ellipses <- df %>%
          dplyr::group_by(.data$cluster) %>%
          dplyr::filter(dplyr::n() >= 3) %>%
          dplyr::ungroup()
        if (nrow(df_ellipses) == 0) {
          return(NULL)
        }
        ggplot2::stat_ellipse(
          data = df_ellipses,
          mapping = ggplot2::aes(x = .data$PC1, y = .data$PC2, group = .data$cluster),
          inherit.aes = FALSE,
          type = "norm",
          level = 0.95,
          linewidth = 0.6,
          linetype = "dashed",
          color = "#4d4d4d",
          show.legend = FALSE
        )
      } else if (boundary_mode == "hull") {
        df_hulls <- df %>%
          dplyr::group_by(.data$cluster) %>%
          dplyr::filter(dplyr::n() >= 3) %>%
          dplyr::group_modify(function(.x, .y) {
            # Compute centroid
            center_x <- mean(.x$PC1)
            center_y <- mean(.x$PC2)

            # Identify hull points
            idx <- grDevices::chull(as.numeric(.x$PC1), as.numeric(.x$PC2))
            hull <- .x[idx, c("PC1", "PC2"), drop = FALSE]

            # Expand hull points radially from centroid
            hull$PC1 <- center_x + (hull$PC1 - center_x) * HULL_EXPANSION
            hull$PC2 <- center_y + (hull$PC2 - center_y) * HULL_EXPANSION

            hull <- dplyr::bind_rows(hull, hull[1, , drop = FALSE])
            # NOTE: Do not return grouping columns (e.g., cluster) from group_modify().
            # dplyr will add them back automatically from the grouping keys.
            tibble::tibble(PC1 = hull$PC1, PC2 = hull$PC2)
          }) %>%
          dplyr::ungroup()
        if (nrow(df_hulls) == 0) {
          return(NULL)
        }
        ggplot2::geom_path(
          data = df_hulls,
          mapping = ggplot2::aes(x = .data$PC1, y = .data$PC2, group = .data$cluster),
          inherit.aes = FALSE,
          linewidth = 0.6,
          linetype = "dashed",
          color = "#4d4d4d",
          show.legend = FALSE
        )
      } else {
        NULL
      }
    } +
    {
      if (show_labels) {
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
        )
      } else {
        NULL
      }
    } +
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
  p <- p + transparent_background_theme()
  ggsave(output_path, p, dpi = 300, width = 10, height = 8, bg = "transparent")
}

write_placeholder_plot <- function(output_path, message_text) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = message_text, size = 4, hjust = 0.5, vjust = 0.5) +
    theme_void()
  p <- p + transparent_background_theme()
  ggsave(output_path, p, dpi = 300, width = 6, height = 4, bg = "transparent")
}
