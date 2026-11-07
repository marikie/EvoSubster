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
