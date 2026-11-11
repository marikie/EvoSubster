# Command-line options, dependencies, and shared plotting constants.

suppressPackageStartupMessages({
  library(argparse)
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(cluster)
  library(scales)
  library(rlang)
  library(ggrepel)
})

EPS <- 1e-12
HULL_EXPANSION <- 1.10
IDENTITY_PREFIX <- "# substitution percent identity:"
DEFAULT_TSV_PATTERN <- "*_ncds.tsv"
CLASSIFICATION_ORDER <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
PHYLUM_COLORS <- c(
  Apicomplexa = "#E69F00",
  Arthropoda = "#0072B2",
  Ascomycota = "#009E73",
  Basidiomycota = "#F0E442",
  Chordata = "#D55E00",
  Cnidaria = "#CC79A7",
  Mucoromycota = "#56B4E9",
  Oomycota = "#000000",
  Unknown = "#999999"
)

invisible(utils::globalVariables(c("mutNum", "totalRootNum", "mutType")))

parse_args <- function() {
  parser <- ArgumentParser(
    description = paste(
      "Traverse dataset directories, collect substitution TSV files from each",
      "dataset's latest completed v2 or latest numeric dated legacy run, and generate 2D PCA and k-means",
      "clustering results. Both log-ratio and observed-frequency variants are produced once. Identity conditions are disabled unless specified; settings and trio decisions are saved as pca_run_settings.json and pca_trio_filter.csv."
    ),
    formatter_class = "argparse.RawDescriptionHelpFormatter",
    epilog = paste(
      "Examples:",
      "  Rscript src/visualize/pca_2d.R results/fungi",
      "  Rscript src/visualize/pca_2d.R results/fungi --output-dir results/fungi/pca_custom",
      "  Rscript src/visualize/pca_2d.R results/fungi --tsv-pattern '*.tsv' --max-clusters 6",
      "  Rscript src/visualize/pca_2d.R results/fungi --cluster-boundary ellipse",
      "  Rscript src/visualize/pca_2d.R results/giw2026 --ingroup-idt-threshold 90 --outgroup-idt-threshold 80 --min-idt-gap-between-out-in 3",
      sep = "\n"
    )
  )
  parser$add_argument(
    "input_dir",
    metavar = "INPUT_DIR",
    help = paste(
      "Root directory containing dataset/trio subdirectories. For each dataset,",
      "the latest completed v2 or numeric dated legacy run is used (for example, results/fungi)."
    )
  )
  parser$add_argument(
    "--tsv-pattern",
    metavar = "GLOB",
    default = DEFAULT_TSV_PATTERN,
    help = sprintf(
      "Supported selection: *_ncds.tsv for noncoding or *.tsv for all genomic single-base channels (default: %s).",
      DEFAULT_TSV_PATTERN
    )
  )
  parser$add_argument(
    "--output-dir",
    metavar = "DIR",
    default = NULL,
    help = "Directory for plots, PCA loadings, clustering metrics, and result CSV files. Default: <INPUT_DIR>/pca."
  )
  parser$add_argument(
    "--max-clusters",
    metavar = "N",
    type = "integer",
    default = 10L,
    help = "Maximum number of k-means clusters to evaluate; the sample count may reduce this limit (default: 10)."
  )
  parser$add_argument(
    "--random-state",
    metavar = "SEED",
    type = "integer",
    default = 42L,
    help = "Random seed used by k-means clustering and reproducible label placement (default: 42)."
  )
  parser$add_argument(
    "--cluster-boundary",
    metavar = "none|ellipse|hull",
    default = "hull",
    help = paste(
      "Boundary overlay for classification-level plots: 'none' draws no boundary,",
      "'ellipse' draws 95%% normal ellipses, and 'hull' draws expanded convex hulls (default: hull)."
    )
  )
  identity_help <- c(
    "idt-threshold" = "Common absolute identity lower bound; individual bounds override it.",
    "ingroup-idt-threshold" = "Ingroup identity (B,C) must strictly exceed this bound.",
    "outgroup-idt-threshold" = "Both outgroup-ingroup identities (A,B and A,C) must strictly exceed this bound.",
    "min-idt-gap-between-out-in" = paste(
      "Both ingroup minus outgroup-ingroup identity gaps must strictly exceed this many percentage points.",
      "Waived for same-genus ingroups and a different-genus outgroup."
    )
  )
  for (option in names(identity_help)) {
    parser$add_argument(
      paste0("--", option), type = "double", default = NULL,
      help = paste(identity_help[[option]], "Optional finite value in [0,100]; disabled when unspecified.")
    )
  }
  args <- parser$parse_args()
  allowed <- c("none", "ellipse", "hull")
  args$cluster_boundary <- tolower(args$cluster_boundary %||% "hull")
  if (!args$cluster_boundary %in% allowed) {
    stop(sprintf("--cluster-boundary must be one of: %s", paste(allowed, collapse = ", ")))
  }
  for (key in c("idt_threshold", "ingroup_idt_threshold", "outgroup_idt_threshold", "min_idt_gap_between_out_in")) {
    value <- args[[key]]
    if (!is.null(value) && (length(value) != 1 || !is.finite(value) || value < 0 || value > 100)) stop(paste(key, "must be finite and in [0,100]"))
  }
  args
}
