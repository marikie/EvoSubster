# Binomial test for each single-base substitution type.
#
# Tests whether each mutType's mutation rate differs significantly from
# the global average rate p = sum(mutNum) / sum(unique oriType totalRootNum).
#
# Usage: Rscript binom_test.R <tsv_file> [label]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript binom_test.R <tsv_file> [label]")
}

tsv_path <- args[1]
label <- if (length(args) >= 2) args[2] else ""

data <- read.delim(tsv_path, header = TRUE, sep = "\t")

# Keep only rows with totalRootNum > 0
data <- data[data$totalRootNum > 0, ]

if (nrow(data) == 0) {
  stop("No rows with totalRootNum > 0")
}

# Derive oriType: "X[Y>Z]W" -> "XYW"
data$oriType <- paste0(substr(data$mutType, 1, 1),
                       substr(data$mutType, 3, 3),
                       substr(data$mutType, nchar(data$mutType), nchar(data$mutType)))

# Global rate: sum(mutNum) / sum(unique oriType totalRootNum)
all_sbst_sum <- sum(data$mutNum)
all_ori_sum <- sum(unique(data[, c("oriType", "totalRootNum")])$totalRootNum)
global_p <- all_sbst_sum / all_ori_sum

# Compute binomial test p-value for each mutType
data$rate <- data$mutNum / data$totalRootNum
data$p_value <- sapply(seq_len(nrow(data)), function(i) {
  binom.test(data$mutNum[i], data$totalRootNum[i], p = global_p)$p.value
})

if (label != "") {
  cat(label, "\n")
}
cat(sprintf("Global rate (p): %.6e\n", global_p))
cat("* p < 0.05, ** p < 0.01\n\n")

# Print results
cat(sprintf("%-12s  %8s  %12s  %12s  %12s  %s\n",
            "mutType", "mutNum", "totalRootNum", "rate", "p_value", "sig"))
cat(paste(rep("-", 68), collapse = ""), "\n")

for (i in seq_len(nrow(data))) {
  sig <- if (data$p_value[i] < 0.01) "**" else if (data$p_value[i] < 0.05) "*" else ""
  cat(sprintf("%-12s  %8d  %12d  %12.6e  %12.4e  %s\n",
              data$mutType[i],
              data$mutNum[i],
              data$totalRootNum[i],
              data$rate[i],
              data$p_value[i],
              sig))
}
