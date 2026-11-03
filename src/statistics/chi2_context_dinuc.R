# Chi-square test for context dependence of dinucleotide substitution rates.
#
# For each original doublet type (e.g., AC, AT, CC, ...),
# tests whether the substitution rate is uniform across all possible
# mutated doublet targets (BC where B != X and C != Y).
# Uses multinomial Monte Carlo (B=10000, seed=1) when expected counts < 5.
# Null probabilities preserve multiplicity when reverse complements are folded.
#
# Usage: Rscript chi2_context_dinuc.R <tsv_file> [label]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript chi2_context_dinuc.R <tsv_file> [label]")
}

tsv_path <- args[1]
label <- if (length(args) >= 2) args[2] else ""

data <- read.delim(tsv_path, header = TRUE, sep = "\t")

required <- c("mutType", "mutNum", "totalRootNum")
if (!all(required %in% names(data)) || anyDuplicated(data$mutType) ||
    any(!is.finite(data$mutNum)) || any(!is.finite(data$totalRootNum)) ||
    any(data$mutNum < 0 | data$mutNum != floor(data$mutNum)) ||
    any(data$totalRootNum < 0 | data$totalRootNum != floor(data$totalRootNum)) ||
    any(data$mutNum > data$totalRootNum)) stop("Invalid DBS count table")
data$orig_doublet <- substr(data$mutType, 1, 2)
reverse_complement <- function(x) {
  paste(chartr("ACGT", "TGCA", rev(strsplit(x, "")[[1]])), collapse = "")
}
excluded <- list(AT = c("GG", "TC", "TG"), CG = c("AA", "AC", "GA"),
                 GC = c("CT", "TG", "TT"), TA = c("AC", "AG", "CC"))
null_weights <- function(root) {
  if (!root %in% c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")) {
    stop("Noncanonical original doublet: ", root)
  }
  targets <- expand.grid(first = strsplit("ACGT", "")[[1]],
                         second = strsplit("ACGT", "")[[1]], stringsAsFactors = FALSE)
  targets <- targets[targets$first != substr(root, 1, 1) &
                       targets$second != substr(root, 2, 2), ]
  folded <- paste0(targets$first, targets$second)
  folded <- vapply(folded, function(x) {
    if (x %in% excluded[[root]]) reverse_complement(x) else x
  }, character(1))
  table(paste0(root, ">", folded)) / 9
}
set.seed(1)
original_doublets <- sort(unique(data$orig_doublet))

if (label != "") {
  cat(label, "\n")
}
cat("* p < 0.05, ** p < 0.01\n")
cat("Null: original nine double substitutions equiprobable; Monte Carlo B=10000, seed=1\n\n")

cat(sprintf("%-8s  %-8s  %12s  %4s  %12s  %s\n",
            "origDblt", "method", "statistic", "df", "p_value", "sig"))
cat(paste(rep("-", 56), collapse = ""), "\n")

for (od in original_doublets) {
  sub <- data[data$orig_doublet == od, ]

  weights <- null_weights(od)
  if (!setequal(sub$mutType, names(weights)) || length(unique(sub$totalRootNum)) != 1 ||
      sum(sub$mutNum) > sub$totalRootNum[1]) stop("Incomplete or inconsistent DBS channels: ", od)
  observed <- sub$mutNum
  probabilities <- as.numeric(weights[sub$mutType])
  total <- sum(observed)
  if (total == 0) {
    cat(sprintf("%-8s  %-8s  %12s  %4s  %12s\n", od, "skip", "NA", "NA", "NA"))
    next
  }
  expected <- total * probabilities
  simulated <- any(expected < 5)
  result <- suppressWarnings(chisq.test(observed, p = probabilities,
                                       simulate.p.value = simulated, B = 10000))
  method <- if (simulated) "multiMC" else "chisq"
  stat <- unname(result$statistic)
  df <- if (simulated) NA else length(observed) - 1
  pval <- result$p.value
  stdres <- (observed - expected) / sqrt(expected * (1 - probabilities))

  sig <- if (pval < 0.01) "**" else if (pval < 0.05) "*" else ""

  stat_str <- if (is.na(stat)) sprintf("%12s", "NA") else sprintf("%12.4f", stat)
  df_str <- if (is.na(df)) sprintf("%4s", "NA") else sprintf("%4d", df)

  cat(sprintf("%-8s  %-8s  %s  %s  %12.4e  %s\n",
              od, method, stat_str, df_str, pval, sig))

  # Print per-context observed, expected counts and standardized residuals
  context_names <- sub$mutType
  cat(sprintf("  %-12s  %8s  %8s  %8s\n", "context", "observed", "expected", "stdres"))
  for (j in seq_along(stdres)) {
    cat(sprintf("  %-12s  %8.0f  %8.2f  %8.2f\n", context_names[j], observed[j], expected[j], stdres[j]))
  }
  cat("\n")
}
