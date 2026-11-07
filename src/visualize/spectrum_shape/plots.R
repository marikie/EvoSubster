shape_save <- function(prefix, draw, width = 14, height = 8) {
  for (ext in c("pdf", "png")) {
    if (ext == "pdf") {
      grDevices::pdf(paste0(prefix, ".pdf"), width = width, height = height)
    } else {
      grDevices::png(paste0(prefix, ".png"), width = width, height = height, units = "in", res = 150)
    }
    tryCatch(draw(), finally = grDevices::dev.off())
  }
}
shape_plots <- function(x, meta, result, embedding, levels, out) {
  annotations <- meta[, ifelse(levels == "species", "taxonomy_species", levels), drop = FALSE]
  names(annotations) <- levels
  annotations$cluster <- as.character(result$membership)
  palettes <- lapply(names(annotations), function(level) shape_palette(annotations[[level]], level))
  names(palettes) <- names(annotations)
  annotation <- ComplexHeatmap::rowAnnotation(df = annotations, col = palettes)
  limit <- max(abs(x))
  if (limit == 0) limit <- 1
  shape_save(file.path(out, "heatmap"), function() {
    ht <- ComplexHeatmap::Heatmap(x,
      name = "log2 ratio", cluster_rows = as.dendrogram(result$tree), row_dend_reorder = FALSE,
      cluster_columns = FALSE, col = circlize::colorRamp2(c(-limit, 0, limit), c("#2166AC", "white", "#B2182B")),
      row_labels = paste(meta$species, meta$display_id), right_annotation = annotation,
      column_names_gp = grid::gpar(fontsize = 6), row_names_gp = grid::gpar(fontsize = 7),
      column_title = "96-channel spectrum similarity (not a phylogeny)"
    )
    ComplexHeatmap::draw(ht)
  }, height = max(8, nrow(x) * 0.2))
  if (!is.null(embedding)) {
    frame <- data.frame(UMAP1 = embedding[, 1], UMAP2 = embedding[, 2], annotations, check.names = FALSE)
    for (level in names(annotations)) {
      frame$group <- frame[[level]]
      p <- ggplot2::ggplot(frame, ggplot2::aes(UMAP1, UMAP2, color = group)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::scale_color_manual(values = palettes[[level]]) +
        ggplot2::labs(color = level) +
        ggplot2::theme_bw()
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
