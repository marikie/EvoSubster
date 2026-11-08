#!/usr/bin/env Rscript
source("src/visualize/spectrum_shape.R")
channels <- shape_channels()
stopifnot(length(unique(channels)) == 96)
df <- data.frame(mutType = channels, mutNum = seq_len(96), totalRootNum = 1000)
x <- shape_vector(df)
stopifnot(isTRUE(all.equal(unname(x), log2(seq_len(96) / mean(seq_len(96))))))
stopifnot(identical(shape_vector(df[96:1, ]), x))
fails <- function(value) stopifnot(inherits(try(value, silent = TRUE), "try-error"))
for (kind in c("zero", "negative", "missing", "duplicate", "extra", "denominator")) {
  z <- df
  if (kind == "zero") z$mutNum[1] <- 0
  if (kind == "negative") z$mutNum[1] <- -1
  if (kind == "missing") z <- z[-1, ]
  if (kind == "duplicate") z$mutType[1] <- z$mutType[2]
  if (kind == "extra") z$mutType[1] <- "bad"
  if (kind == "denominator") z$totalRootNum[1] <- 0
  fails(shape_vector(z))
}
x <- rbind(a = x, b = x + 0.01, c = -x, d = -x + 0.01)
r <- shape_cluster(x)
stopifnot(r$k == 2, r$membership[1] == r$membership[2], r$membership[1] != r$membership[3])
stopifnot(identical(unname(r$representatives), c("a", "c")))
stopifnot(shape_cluster(x[1:2, ])$k == 1)
stopifnot(shape_cluster(x * 0)$k == 1)
stopifnot(shape_cluster(x, clusters = "3")$k == 3)
fails(shape_cluster(x, clusters = "5"))
fails(shape_cluster(x * 0, method = "correlation"))
stopifnot(is.null(shape_umap(x[1:3, ], "ward", 42)))
cat("Spectrum shape core tests passed\n")
stopifnot(shape_cluster(x, method = "correlation")$k >= 1)

correlation <- shape_cluster(x, method = "correlation", clusters = "2")
stopifnot(isTRUE(all.equal(as.matrix(correlation$distance), 1 - cor(t(x)), tolerance = 1e-12)))
stopifnot(correlation$membership[1] == correlation$membership[2], correlation$membership[1] != correlation$membership[3])

annotations <- data.frame(phylum = c("Cnidaria", "Unknown", "Mucoromycota"), class = c("Unknown", "Same", "Other"), order = c("Same", "Unknown", "Other"), cluster = c("1", "2", "3"))
colors <- shape_annotation_colors(annotations)
stopifnot(!anyDuplicated(colors$color), nrow(colors) == 12)
stopifnot(identical(colors, shape_annotation_colors(annotations[3:1, ])))
stopifnot(colors$color[colors$level == "phylum" & colors$label == "Cnidaria"] == "#CC79A7")
stopifnot(length(unique(colors$color[colors$label == "Unknown"])) == 3)

stopifnot(!anyDuplicated(colors$identifier))
stopifnot(identical(colors$identifier[colors$level == "class"], c("C1", "C2", "C3")))
stopifnot(identical(colors$identifier[colors$level == "cluster"], c("K1", "K2", "K3")))
extra <- shape_annotation_colors(data.frame(genus = c("Z", "A"), cluster = c("10", "2")))
stopifnot(identical(extra$identifier[extra$level == "genus"], c("genus1", "genus2")))
stopifnot(identical(extra$identifier[extra$level == "cluster"], c("K10", "K2")))
runs <- shape_annotation_runs(c("a", "b", "a", "b", "a"), c(5, 1, 2, 4, 3))
stopifnot(identical(runs$label, c("a", "b", "a")), identical(runs$start, c(1L, 3L, 5L)))
stopifnot(identical(runs$end, c(2L, 4L, 5L)), identical(runs$center, c(1.5, 3.5, 5)))
stopifnot(isTRUE(all.equal(runs$y, c(0.8, 0.4, 0.1))))
stopifnot(identical(shape_contrast_text(c("#000000", "#FFFFFF")), c("white", "black")))
coord <- shape_color_coordinates(c("#999999", "#FF0000", "#00FF00", "#0000FF"))
stopifnot(is.na(coord$hue[1]), all(is.finite(coord$hue[-1])))
stopifnot(shape_score_better(c(-1, 1, 1), c(-2, 100, 100)))
stopifnot(shape_score_better(c(0, 11, 1), c(0, 10, 100)))
stopifnot(shape_score_better(c(0, 10, 11), c(0, 10, 10)))
stopifnot(!shape_score_better(c(0, 10, 10), c(0, 10, 10)))
# An actual class/order adjacency overrides their initially similar purple hues.
close <- data.frame(phylum = c("Basidiomycota", "Chordata"),
  class = c("Agaricomycetes", "Mammalia"), order = c("Agaricales", "Primates"), cluster = c("1", "2"))
optimized <- shape_annotation_colors(close)
target <- optimized$color[match(c("Agaricomycetes", "Agaricales"), optimized$label)]
hues <- shape_color_coordinates(target)$hue
gap <- abs(diff(hues))
stopifnot(min(gap, 360 - gap) >= 60)
stopifnot(optimized$color[optimized$label == "Basidiomycota"] == "#F0E442")
stopifnot(identical(optimized, shape_annotation_colors(close[2:1, ])))
# Adjacency is defined by displayed columns, not all pairs or row frequency.
edges <- shape_annotation_edges(close, optimized)
stopifnot(nrow(edges) == 6, identical(edges, shape_annotation_edges(close[c(2, 1, 1), ], optimized)))
cat("Spectrum annotation adjacency, identifiers and run-placement tests passed\n")

# Brackets use displayed contiguous runs, including separated repeats.
layout <- shape_cluster_bracket_layout(c("1", "10000", "1", "2"), c(1, 3, 2, 4), 40, c(12, 20, 12))
stopifnot(identical(layout$runs$label, c("1", "10000", "2")))
stopifnot(identical(layout$runs$top, c(1, 0.5, 0.25)))
stopifnot(identical(layout$runs$bottom, c(0.5, 0.25, 0)))
stopifnot(identical(layout$runs$text, c("Cluster 1", "Cluster 10000", "Cluster 2")))
stopifnot(all(layout$runs$fontsize >= 6 & layout$runs$fontsize <= 8), !all(layout$runs$fits))
grown <- shape_cluster_bracket_layout(c("1", "10000", "1", "2"), c(1, 3, 2, 4), layout$required_body_mm, c(12, 20, 12))
stopifnot(all(grown$runs$fits), grown$runs$fontsize[2] == 6)
singleton <- shape_cluster_bracket_layout("123456789", 1L, 8, 30)
stopifnot(!singleton$runs$fits, singleton$required_body_mm == 24.5)
stopifnot(shape_cluster_bracket_layout("123456789", 1L, singleton$required_body_mm, 30)$runs$fits)
repeated <- shape_cluster_bracket_layout(c("1", "2", "1"), 1:3, 60, rep(12, 3))
stopifnot(identical(repeated$runs$label, c("1", "2", "1")), all(repeated$runs$fits))
cat("Cluster bracket geometry and short-run fitting tests passed\n")
