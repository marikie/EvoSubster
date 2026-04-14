# Chi-square test for context dependence of dinucleotide substitution rates.
#
# For each original doublet type (e.g., AC, AT, CC, ...),
# tests whether the substitution rate is uniform across all possible
# mutated doublet targets (BC where B != X and C != Y).
# Falls back to Fisher's exact test (Monte Carlo) when expected counts < 5.
#
# Usage: Rscript chi2_context_dinuc.R <tsv_file> [label]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript chi2_context_dinuc.R <tsv_file> [label]")
}

tsv_path <- args[1]
label <- if (length(args) >= 2) args[2] else ""

data <- read.delim(tsv_path, header = TRUE, sep = "\t")

# Extract original doublet: "XY>BC" -> "XY"
data$orig_doublet <- substr(data$mutType, 1, 2)

# Keep only rows with totalRootNum > 0
data <- data[data$totalRootNum > 0, ]

if (nrow(data) == 0) {
  stop("No rows with totalRootNum > 0")
}

original_doublets <- sort(unique(data$orig_doublet))

if (label != "") {
  cat(label, "\n")
}
cat("* p < 0.05, ** p < 0.01\n\n")

cat(sprintf("%-8s  %-8s  %12s  %4s  %12s  %s\n",
            "origDblt", "method", "statistic", "df", "p_value", "sig"))
cat(paste(rep("-", 56), collapse = ""), "\n")

for (od in original_doublets) {
  sub <- data[data$orig_doublet == od, ]

  # Build 2 x N contingency table: row 1 = mutNum, row 2 = totalRootNum - mutNum
  mut_counts <- sub$mutNum
  non_counts <- sub$totalRootNum - sub$mutNum
  mat <- rbind(mut_counts, non_counts)

  # Remove columns where both rows are zero
  keep <- colSums(mat) > 0
  mat <- mat[, keep, drop = FALSE]

  if (ncol(mat) < 2) {
    cat(sprintf("%-8s  %-8s  %12s  %4s  %12s  %s\n",
                od, "skip", "NA", "NA", "NA", ""))
    next
  }

  # Try chi-square test; check for low expected counts
  exp_counts <- suppressWarnings(chisq.test(mat))$expected
  has_low_expected <- any(exp_counts < 5)

  if (!has_low_expected) {
    result <- chisq.test(mat)
    method <- "chisq"
    stat <- result$statistic
    df <- result$parameter
    pval <- result$p.value
    stdres <- result$stdres[1, ]  # standardized residuals for mut row
    expected <- result$expected[1, ]  # expected counts for mut row
  } else {
    result <- fisher.test(mat, simulate.p.value = TRUE, B = 10000)
    method <- "fisher"
    stat <- NA
    df <- NA
    pval <- result$p.value
    # Compute stdres and expected manually for fisher fallback
    chi_result <- suppressWarnings(chisq.test(mat))
    stdres <- chi_result$stdres[1, ]
    expected <- chi_result$expected[1, ]
  }

  sig <- if (pval < 0.01) "**" else if (pval < 0.05) "*" else ""

  stat_str <- if (is.na(stat)) sprintf("%12s", "NA") else sprintf("%12.4f", stat)
  df_str <- if (is.na(df)) sprintf("%4s", "NA") else sprintf("%4d", df)

  cat(sprintf("%-8s  %-8s  %s  %s  %12.4e  %s\n",
              od, method, stat_str, df_str, pval, sig))

  # Print per-context observed, expected counts and standardized residuals
  context_names <- sub$mutType[keep]
  observed <- mat[1, ]
  cat(sprintf("  %-12s  %8s  %8s  %8s\n", "context", "observed", "expected", "stdres"))
  for (j in seq_along(stdres)) {
    cat(sprintf("  %-12s  %8.0f  %8.2f  %8.2f\n", context_names[j], observed[j], expected[j], stdres[j]))
  }
  cat("\n")
}
