#!/usr/bin/env Rscript
# Filter chi2_context output: show only sections where the specified context
# has the highest standardized residual for the given mutation type.
#
# Usage: Rscript filter_top_context.R <chi2_output_file> <mut_type> <context>
# Example:
#   Rscript filter_top_context.R fungi_chi2_context.txt "T>C" "T[T>C]A"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript filter_top_context.R <chi2_output_file> <mut_type> <context>")
}

file_path  <- args[1]
target_mut <- args[2]
target_ctx <- args[3]

lines <- readLines(file_path, warn = FALSE)

# Parse lines into records: (section, mutType, p_value, context, stdres)
records         <- list()
current_section <- NA_character_
current_mut     <- NA_character_
current_pval    <- NA_real_

for (line in lines) {
  # Section header: === name ===
  if (grepl("^=== .+ ===$", line)) {
    current_section <- sub("^=== (.+) ===$", "\\1", line)
    current_mut     <- NA_character_
    current_pval    <- NA_real_
    next
  }

  if (is.na(current_section)) next

  # Mutation type summary line, e.g.: "T>C  chisq  455.01  15  1.80e-87  **"
  if (grepl("^(C>A|C>G|C>T|T>A|T>C|T>G)\\s", line)) {
    tokens       <- strsplit(trimws(line), "\\s+")[[1]]
    current_mut  <- tokens[1]
    current_pval <- suppressWarnings(as.numeric(tokens[5]))
    next
  }

  # Context/stdres line, e.g.: "  A[C>A]A    9.10"
  if (!is.na(current_mut) && grepl("^  [ACGT]\\[", line)) {
    tokens <- strsplit(trimws(line), "\\s+")[[1]]
    records[[length(records) + 1]] <- list(
      section = current_section,
      mutType = current_mut,
      p_value = current_pval,
      context = tokens[1],
      stdres  = suppressWarnings(as.numeric(tokens[2]))
    )
  }
}

if (length(records) == 0) {
  cat("No data parsed from file.\n")
  quit(status = 0)
}

df <- do.call(rbind, lapply(records, as.data.frame, stringsAsFactors = FALSE))
df$stdres  <- as.numeric(df$stdres)
df$p_value <- as.numeric(df$p_value)

# Restrict to target mutation type, drop NA stdres (skip entries)
sub <- df[df$mutType == target_mut & !is.na(df$stdres), ]

if (nrow(sub) == 0) {
  cat(sprintf("No data found for mutation type '%s'.\n", target_mut))
  quit(status = 0)
}

# Per-section maximum stdres
max_by_section <- tapply(sub$stdres, sub$section, max)

# Keep sections where target_ctx achieves the maximum
hits <- sub[sub$context == target_ctx &
              sub$stdres == max_by_section[sub$section], ]

cat(sprintf("Filter: %s, context = %s (highest stdres in category)\n\n",
            target_mut, target_ctx))

if (nrow(hits) == 0) {
  cat("No matching sections found.\n")
  quit(status = 0)
}

hits <- hits[order(hits$stdres, decreasing = TRUE), ]

sig_col <- ifelse(is.na(hits$p_value), "",
           ifelse(hits$p_value < 0.01, "**",
           ifelse(hits$p_value < 0.05, "*", "")))

cat(sprintf("%-52s  %12s  %-2s  %8s\n", "section", "p_value", "sig", "stdres"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(hits))) {
  cat(sprintf("%-52s  %12.4e  %-2s  %8.2f\n",
              hits$section[i], hits$p_value[i], sig_col[i], hits$stdres[i]))
}
cat("\n")
