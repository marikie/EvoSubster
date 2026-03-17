#! /usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1 || length(args) > 2) {
  stop("Usage: dinucleotide-plot.R <infile> [out_dir]")
}

inFile <- args[1]
out_dir <- if (length(args) >= 2) args[2] else dirname(inFile)
outFile <- file.path(out_dir, paste(basename(inFile), "pdf", sep="."))

xlab <- "Substitutions / original dinucleotides (%)"

d <- read.table(inFile, header=TRUE)
barLabels <- d[[1]]
percentages <- 100 * d[[2]] / d[[3]]

pdf(outFile, pointsize=8)

par(mar=c(2.6, 2.1, 0.1, 1))  # margin widths
par(mgp=c(1.5, 0.4, 0))       # axis label positions
par(tcl=-0.3)                 # length of axis tick marks
par(yaxs="i")
par(las=1)

# colors:
col <- c(rep(2,9), rep(3,6), rep(4,9), rep(5,6), rep(6,9),
         rep(7,6), rep(8,6), rep(2,9), rep(3,9), rep(4,9))

pos <- barplot(percentages, space=0, horiz=TRUE, col=col, xlab=xlab)

axis(2, pos, barLabels, tick=FALSE, line=1.5, hadj=0, cex.axis=0.8)
