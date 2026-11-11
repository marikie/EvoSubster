shape_save <- function(prefix, draw, width = 14, height = 8) {
  for (ext in c("pdf", "png")) {
    if (ext == "pdf") {
      grDevices::pdf(paste0(prefix, ".pdf"), width = width, height = height, bg = "transparent")
    } else {
      grDevices::png(paste0(prefix, ".png"), width = width, height = height, units = "in", res = 150, bg = "transparent")
    }
    tryCatch(draw(), finally = grDevices::dev.off())
  }
}
shape_plots <- function(x, meta, result, embedding, levels, out) {
  annotations <- meta[, ifelse(levels == "species", "taxonomy_species", levels), drop = FALSE]
  names(annotations) <- levels
  annotations$cluster <- as.character(result$membership)
  colors <- shape_annotation_colors(annotations)
  shape_csv(colors, file.path(out, "annotation_colors.csv"))
  palettes <- lapply(names(annotations), function(level) {
    rows <- colors[colors$level == level, ]
    setNames(rows$color, rows$label)
  })
  names(palettes) <- names(annotations)
  identifiers <- lapply(names(annotations), function(level) {
    rows <- colors[colors$level == level, ]
    setNames(rows$identifier, rows$label)
  })
  names(identifiers) <- names(annotations)
  limit <- max(abs(x))
  if (limit == 0) limit <- 1
  bracket_measure <- new.env(parent = emptyenv())
  draw_heatmap <- function() {
    legends <- lapply(levels, function(level) {
      labels <- names(palettes[[level]])
      ComplexHeatmap::Legend(title = level, at = labels,
        labels = paste0(identifiers[[level]][labels], ": ", labels),
        legend_gp = grid::gpar(fill = palettes[[level]][labels]))
    })
    strips <- lapply(levels, function(level) {
      values <- as.character(annotations[[level]])
      palette <- palettes[[level]]
      ids <- identifiers[[level]]
      ComplexHeatmap::AnnotationFunction(which = "row", n = length(values), width = grid::unit(8, "mm"),
        var_import = list(values = values, palette = palette, ids = ids),
        fun = function(index, k = 1, n = 1) {
          count <- length(index)
          grid::grid.rect(x = 0.5, y = 1 - (seq_len(count) - 0.5) / count,
            width = 1, height = 1 / count, gp = grid::gpar(fill = palette[values[index]], col = NA))
          runs <- shape_annotation_runs(values, index)
          labels <- unname(ids[runs$label])
          # Fit longer prefixes for additional ranks inside the same 8 mm strip.
          grid::grid.text(labels, x = 0.5, y = runs$y,
            gp = grid::gpar(col = shape_contrast_text(palette[runs$label]),
              fontsize = pmin(8, 32 / nchar(labels)), fontface = "bold"))
        })
    })
    names(strips) <- levels
    row_labels <- paste(meta$species, meta$display_id)
    # Explicit text annotation keeps brackets after full species labels.
    label_gp <- grid::gpar(fontsize = 7)
    strips$sample_labels <- ComplexHeatmap::anno_text(row_labels, which = "row",
      location = grid::unit(0, "mm"), just = "left", gp = label_gp,
      width = ComplexHeatmap::max_text_width(row_labels, gp = label_gp))
    cluster_values <- as.character(result$membership)
    strips$cluster_brackets <- ComplexHeatmap::AnnotationFunction(which = "row", n = nrow(x),
      width = grid::unit(8, "mm"), var_import = list(cluster_values = cluster_values, bracket_measure = bracket_measure),
      fun = function(index, k = 1, n = 1) {
        runs <- shape_annotation_runs(cluster_values, index)
        labels <- paste("Cluster", runs$label)
        text_width <- vapply(labels, function(label) {
          grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = 8))), "mm", valueOnly = TRUE)
        }, numeric(1))
        body <- grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
        layout <- shape_cluster_bracket_layout(cluster_values, index, body, text_width)
        bracket_measure$body_mm <- body
        bracket_measure$required_body_mm <- layout$required_body_mm
        runs <- layout$runs
        # ] brackets open toward the species labels; all share the same x.
        grid::grid.segments(x0 = grid::unit(1.5, "mm"), x1 = grid::unit(1.5, "mm"),
          y0 = runs$bottom, y1 = runs$top, gp = grid::gpar(col = "grey55", lwd = 0.8))
        grid::grid.segments(x0 = grid::unit(0, "mm"), x1 = grid::unit(1.5, "mm"),
          y0 = c(runs$bottom, runs$top), y1 = c(runs$bottom, runs$top),
          gp = grid::gpar(col = "grey55", lwd = 0.8))
        grid::grid.text(runs$text, x = grid::unit(5, "mm"), y = runs$y, rot = 90,
          gp = grid::gpar(col = "black", fontsize = runs$fontsize))
      })
    annotation <- do.call(ComplexHeatmap::rowAnnotation, c(strips,
      list(gap = grid::unit(2, "mm"), show_annotation_name = c(rep(TRUE, length(levels)), FALSE, FALSE))))
    ht <- ComplexHeatmap::Heatmap(x,
      name = "log2 ratio", cluster_rows = as.dendrogram(result$tree), row_dend_reorder = FALSE,
      cluster_columns = FALSE, col = circlize::colorRamp2(c(-limit, 0, limit), c("#2166AC", "white", "#B2182B")),
      show_row_names = FALSE, right_annotation = annotation,
      column_names_gp = grid::gpar(fontsize = 6, fontfamily = "mono"), row_names_gp = grid::gpar(fontsize = 7),
      column_title = "96-channel spectrum similarity (not a phylogeny)"
    )
    ComplexHeatmap::draw(ht, annotation_legend_list = legends, background = "transparent")
  }
  heatmap_height <- max(8, nrow(x) * 0.2)
  # Measure real layout without creating a file, then grow for short cluster runs.
  grDevices::pdf(NULL, width = 15, height = heatmap_height, bg = "transparent")
  tryCatch(draw_heatmap(), finally = grDevices::dev.off())
  heatmap_height <- heatmap_height + max(0, bracket_measure$required_body_mm - bracket_measure$body_mm) / 25.4
  shape_save(file.path(out, "heatmap"), draw_heatmap, width = 15, height = heatmap_height)
  if (!is.null(embedding)) {
    frame <- data.frame(UMAP1 = embedding[, 1], UMAP2 = embedding[, 2], annotations, check.names = FALSE)
    for (level in names(annotations)) {
      frame$group <- frame[[level]]
      p <- ggplot2::ggplot(frame, ggplot2::aes(UMAP1, UMAP2, color = group)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::scale_color_manual(values = palettes[[level]],
          breaks = names(palettes[[level]]),
          labels = paste0(identifiers[[level]], ": ", names(palettes[[level]]))) +
        ggplot2::labs(color = level) +
        ggplot2::theme_bw() + transparent_background_theme()
      shape_save(file.path(out, paste0("umap_", level)), function() print(p), 10, 7)
    }
  }
  summaries <- list()
  representatives <- list()
  for (g in sort(unique(result$membership))) {
    member <- x[result$membership == g, , drop = FALSE]
    representative <- result$representatives[g]
    q <- apply(member, 2, quantile, probs = c(0.25, 0.5, 0.75))
    summaries[[g]] <- data.frame(cluster = g, channel = colnames(x), q25 = q[1, ], median = q[2, ], q75 = q[3, ], n = nrow(member))
    representatives[[g]] <- data.frame(cluster = g, sample_id = representative, channel = colnames(x), log2_ratio = x[representative, ])
  }
  summary <- do.call(rbind, summaries)
  reps <- do.call(rbind, representatives)
  shape_csv(summary, file.path(out, "cluster_summaries.csv"))
  shape_csv(reps, file.path(out, "representative_spectra.csv"))
  range <- range(c(x, summary$q25, summary$q75, 0))
  if (diff(range) == 0) range <- c(-1, 1)
  for (g in sort(unique(result$membership))) {
    s <- summary[summary$cluster == g, ]
    r <- reps[reps$cluster == g, ]
    members <- meta[result$membership == g, , drop = FALSE]
    composition <- paste(vapply(levels, function(level) paste0(level, ": ", paste(names(table(members[[if (level == "species") "taxonomy_species" else level]])), table(members[[if (level == "species") "taxonomy_species" else level]]), collapse = "; ")), character(1)), collapse = "\n")
    shape_save(file.path(out, paste0("spectrum_cluster_", g)), function() {
      par(mfrow = c(2, 1), mar = c(6, 4, 3, 1), oma = c(3, 0, 0, 0))
      bars <- barplot(r$log2_ratio,
        col = rep(c("#55AACC", "#222222", "#CC4444", "#AAAAAA", "#66AA55", "#DD88AA"), each = 16), ylim = range, border = NA, names.arg = r$channel, las = 2, cex.names = 0.45,
        main = paste0("Cluster ", g, " (n=", nrow(members), "): ", members$species[match(r$sample_id[1], members$sample_id)]), ylab = "log2 ratio"
      )
      abline(h = 0, col = "grey")
      plot(seq_len(96), s$median, type = "n", ylim = range, xaxt = "n", xlab = "", ylab = "log2 ratio", main = "Median and interquartile range")
      segments(seq_len(96), s$q25, seq_len(96), s$q75, col = "grey50")
      points(seq_len(96), s$median, pch = 16, cex = 0.5)
      axis(1, at = seq_len(96), labels = s$channel, las = 2, cex.axis = 0.45)
      abline(h = 0, col = "grey")
      mtext(composition, side = 1, outer = TRUE, cex = 0.65)
    }, 16, 10)
  }
}
