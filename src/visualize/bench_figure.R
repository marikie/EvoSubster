#!/usr/bin/env Rscript

# bench_figure.R
# Generate supplementary figure and table for computational performance of
# the alignment + substitution-counting pipeline.
#
# Usage:
#   Rscript bench_figure.R <bench_summary.tsv> <bench_metadata.tsv> <out_dir>
#
# Outputs:
#   <out_dir>/bench_figure.pdf  — two-panel scatter plot
#   <out_dir>/bench_table.tsv   — formatted summary table

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)
library(RColorBrewer)

# ---- parameters -------------------------------------------------------
POINT_SIZE   <- 2.5
POINT_ALPHA  <- 0.85
LABEL_SIZE   <- 2.0    # ggrepel label size (pt)
LINE_COLOR   <- "grey40"
LINE_WIDTH   <- 0.5
LINEAGE_COLORS <- c(cnidaria = "#2166AC", fungi = "#D6604D")
FIGURE_WIDTH  <- 22   # cm
FIGURE_HEIGHT <- 40   # cm (two panels stacked vertically)
# -----------------------------------------------------------------------

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: bench_figure.R <bench_summary.tsv> <bench_metadata.tsv> <out_dir>")
  }
  list(bench_summary = args[1], bench_metadata = args[2], out_dir = args[3])
}

# Parse wall-clock string (MM:SS.ss or H:MM:SS) to decimal minutes.
parse_wall_clock_min <- function(x) {
  sapply(x, function(s) {
    parts <- strsplit(s, ":")[[1]]
    if (length(parts) == 2) {
      as.numeric(parts[1]) + as.numeric(parts[2]) / 60
    } else {
      as.numeric(parts[1]) * 60 + as.numeric(parts[2]) + as.numeric(parts[3]) / 60
    }
  }, USE.NAMES = FALSE)
}

main <- function() {
  args <- parse_args()
  dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- load data -------------------------------------------------------
  bench <- read.table(args$bench_summary, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "")
  meta  <- read.table(args$bench_metadata, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "")

  # ---- join and compute derived columns --------------------------------
  df <- inner_join(bench, meta, by = "trio_label") %>%
    mutate(
      runtime_min  = parse_wall_clock_min(wall_clock),
      total_mbp    = (as.numeric(size_org1_bp) + as.numeric(size_org2_bp) + as.numeric(size_org3_bp)) / 1e6,
      min_identity = pmin(idt_12, idt_13, idt_23),
      peak_ram_gb  = max_rss_kb / 1e6,
      # Short label for labeling extreme points
      short_label  = trio_label
    ) %>%
    # Use lineage from bench_summary (column lineage.x after join)
    rename(lineage = lineage.x) %>%
    select(-lineage.y)

  # ---- write table -----------------------------------------------------
  table_out <- df %>%
    mutate(
      org1_mbp  = round(size_org1_bp / 1e6, 1),
      org2_mbp  = round(size_org2_bp / 1e6, 1),
      org3_mbp  = round(size_org3_bp / 1e6, 1),
      total_mbp = round(total_mbp, 1),
      min_idt   = round(min_identity, 2),
      ram_gb    = round(peak_ram_gb, 2)
    ) %>%
    arrange(desc(runtime_min)) %>%
    select(
      trio_label, lineage,
      org1_mbp, org2_mbp, org3_mbp, total_mbp,
      min_idt, wall_clock, cpu_percent, ram_gb
    )

  colnames(table_out) <- c(
    "trio", "lineage",
    "org1_size_Mbp", "org2_size_Mbp", "org3_size_Mbp", "total_size_Mbp",
    "min_pairwise_identity_pct", "wall_clock", "cpu_pct", "peak_RAM_GB"
  )

  write.table(table_out,
              file = file.path(args$out_dir, "bench_table.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  message("Wrote bench_table.tsv")

  repel_opts <- list(
    size               = LABEL_SIZE,
    show.legend        = FALSE,
    max.overlaps       = Inf,
    box.padding        = 0.5,
    point.padding      = 0.3,
    force              = 8,
    force_pull         = 0.3,
    min.segment.length = 0.2,
    max.iter           = 50000,
    seed               = 42,
    fill               = "white",
    label.size         = NA,
    label.padding      = unit(0.12, "lines"),
    alpha              = 0.9
  )

  # ---- plot: Panel A — genome size vs runtime --------------------------
  pA <- ggplot(df, aes(x = total_mbp, y = runtime_min, color = lineage)) +
    geom_smooth(method = "lm", se = FALSE, color = LINE_COLOR,
                linewidth = LINE_WIDTH, linetype = "dashed") +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    do.call(geom_label_repel, c(list(aes(label = short_label)), repel_opts)) +
    scale_x_log10(limits = c(30, NA), labels = scales::label_comma()) +
    scale_y_log10(limits = c(2.5, NA), labels = scales::label_comma()) +
    scale_color_manual(values = LINEAGE_COLORS) +
    coord_cartesian(clip = "off") +
    labs(
      x     = "Total assembly size (Mbp, log scale)",
      y     = "Wall-clock time (min, log scale)",
      color = "Lineage",
      tag   = "A"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.tag        = element_text(face = "bold"),
      plot.margin     = margin(5, 10, 10, 10)
    )

  # ---- plot: Panel B — min pairwise identity vs runtime ----------------
  pB <- ggplot(df, aes(x = min_identity, y = runtime_min, color = lineage)) +
    geom_smooth(method = "lm", se = FALSE, color = LINE_COLOR,
                linewidth = LINE_WIDTH, linetype = "dashed") +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    do.call(geom_label_repel, c(list(aes(label = short_label)), repel_opts)) +
    scale_x_continuous(limits = c(78, NA)) +
    scale_y_log10(limits = c(2.5, NA), labels = scales::label_comma()) +
    scale_color_manual(values = LINEAGE_COLORS) +
    coord_cartesian(clip = "off") +
    labs(
      x     = "Minimum pairwise sequence identity (%)",
      y     = "Wall-clock time (min, log scale)",
      color = "Lineage",
      tag   = "B"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.tag        = element_text(face = "bold"),
      plot.margin     = margin(5, 10, 5, 10)
    )

  # ---- combine and save ------------------------------------------------
  combined <- pA / pB + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  out_pdf <- file.path(args$out_dir, "bench_figure.pdf")
  pdf(out_pdf, width = FIGURE_WIDTH / 2.54, height = FIGURE_HEIGHT / 2.54)
  print(combined)
  dev.off()
  message("Wrote bench_figure.pdf")
}

main()
