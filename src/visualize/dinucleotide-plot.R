#! /usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1 || length(args) > 3) {
  stop("Usage: dinucleotide-plot.R <infile> [out_dir] [species_full_name]")
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "plot_species_label.R"))

inFile <- args[1]
out_dir <- if (length(args) >= 2) args[2] else dirname(inFile)
outFile <- file.path(out_dir, paste(basename(inFile), "pdf", sep="."))
species_name <- plot_species_name(inFile, if (length(args) >= 3) args[3] else NULL)

xlab <- "Substitutions / original dinucleotides (%)"

d <- read.table(inFile, header=TRUE)
barLabels <- d[[1]]
percentages <- 100 * d[[2]] / d[[3]]

# Mark group maxima that contain at least one transversion.
original <- sub(">.*$", "", barLabels)
mutated <- sub("^.*>", "", barLabels)
transitions <- c("CT", "TC", "GA", "AG")
has_transversion <- rep(FALSE, length(barLabels))
for (i in seq_len(2)) {
  from <- substr(original, i, i)
  to <- substr(mutated, i, i)
  has_transversion <- has_transversion |
    (from != to & !(paste0(from, to) %in% transitions))
}
group_max <- ave(percentages, original, FUN = max)
marked <- which(is.finite(percentages) & percentages == group_max & has_transversion)

pdf(outFile, pointsize=8)

par(mar=c(2.6, 4.5, 2, 1))    # margin widths, including species label
par(mgp=c(1.5, 0.4, 0))       # axis label positions
par(tcl=-0.3)                 # length of axis tick marks
par(yaxs="i")
par(las=1)

# colors:
col <- c(rep(2,9), rep(3,6), rep(4,9), rep(5,6), rep(6,9),
         rep(7,6), rep(8,6), rep(2,9), rep(3,9), rep(4,9))

pos <- barplot(percentages, space=0, horiz=TRUE, col=col, xlab=xlab,
               xlim=c(0, max(percentages, na.rm=TRUE) * 1.05))

# Draw five-point stars about 5 points wide, with a 2-point gap from each bar.
# Convert physical sizes to plot coordinates so the stars keep their shape.
star_angles <- pi / 2 + (0:9) * pi / 5
star_radii <- rep(c(2.5, 1), 5) / 72
star_x <- xinch(star_radii * cos(star_angles))
star_y <- yinch(star_radii * sin(star_angles))
star_offset <- xinch(2 / 72) - min(star_x)
for (i in marked) {
  polygon(percentages[i] + star_offset + star_x, pos[i] + star_y,
          col="black", border=NA)
}

axis(2, pos, barLabels, tick=FALSE, line=0, hadj=1, cex.axis=0.8, family="mono")
add_species_label(species_name, cex = 0.9)
dev.off()
