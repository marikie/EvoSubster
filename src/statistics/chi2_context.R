# Chi-square test for context dependence of substitution rates.
#
# For each of 6 central mutation types (C>A, C>G, C>T, T>A, T>C, T>G),
# tests whether the substitution rate is uniform across 16 flanking contexts.
# Falls back to Fisher's exact test (Monte Carlo) when expected counts < 5.
#
# Usage: Rscript chi2_context.R <tsv_file> [label]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript chi2_context.R <tsv_file> [label]")
}

tsv_path <- args[1]
label <- if (length(args) >= 2) args[2] else ""

data <- read.delim(tsv_path, header = TRUE, sep = "\t")

# Extract central mutation type: "X[Y>Z]W" -> "Y>Z"
data$central <- substr(data$mutType, 3, 5)

# Keep only rows with totalRootNum > 0
data <- data[data$totalRootNum > 0, ]

if (nrow(data) == 0) {
  stop("No rows with totalRootNum > 0")
}

central_types <- sort(unique(data$central))

if (label != "") {
  cat(label, "\n")
}
cat("* p < 0.05, ** p < 0.01\n\n")

cat(sprintf("%-8s  %-8s  %12s  %4s  %12s  %s\n",
            "mutType", "method", "statistic", "df", "p_value", "sig"))
cat(paste(rep("-", 56), collapse = ""), "\n")

for (ct in central_types) {
  sub <- data[data$central == ct, ]

  # Build 2 x N contingency table: row 1 = mutNum, row 2 = totalRootNum - mutNum
  mut_counts <- sub$mutNum
  non_counts <- sub$totalRootNum - sub$mutNum
  mat <- rbind(mut_counts, non_counts)

  # Remove columns where both rows are zero
  keep <- colSums(mat) > 0
  mat <- mat[, keep, drop = FALSE]

  if (ncol(mat) < 2) {
    cat(sprintf("%-8s  %-8s  %12s  %4s  %12s  %s\n",
                ct, "skip", "NA", "NA", "NA", ""))
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
              ct, method, stat_str, df_str, pval, sig))

  # Print per-context expected counts and standardized residuals
  context_names <- sub$mutType[keep]
  cat(sprintf("  %-12s  %8s  %8s\n", "context", "expected", "stdres"))
  for (j in seq_along(stdres)) {
    cat(sprintf("  %-12s  %8.2f  %8.2f\n", context_names[j], expected[j], stdres[j]))
  }
  cat("\n")
}
