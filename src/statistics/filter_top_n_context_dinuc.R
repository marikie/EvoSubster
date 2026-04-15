#!/usr/bin/env Rscript
# Filter chi2_context_dinuc output: show only sections where ALL specified
# contexts rank within the top N standardized residuals for the given
# original doublet type, where N = number of specified contexts.
#
# Input file must contain one or more sections delimited by === name === headers,
# each section being the output of chi2_context_dinuc.R.
#
# Usage: Rscript filter_top_n_context_dinuc.R <chi2_output_file> <orig_doublet> <ctx1> <ctx2> [ctx3 ...]
# Example:
#   Rscript filter_top_n_context_dinuc.R glom_chi2_dinuc.txt "AC" "AC>GT" "AC>TT"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript filter_top_n_context_dinuc.R <chi2_output_file> <orig_doublet> <ctx1> [ctx2 ...]")
}

file_path   <- args[1]
target_mut  <- args[2]
target_ctxs <- args[3:length(args)]
K           <- length(target_ctxs)

lines <- readLines(file_path, warn = FALSE)

# Parse lines into records: (section, origDblt, p_value, context, stdres)
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

  # Doublet summary line, e.g.: "AC  chisq  1566.9453  8  0.0000e+00  **"
  if (grepl("^[ACGT]{2}\\s", line)) {
    tokens       <- strsplit(trimws(line), "\\s+")[[1]]
    current_mut  <- tokens[1]
    current_pval <- suppressWarnings(as.numeric(tokens[5]))
    next
  }

  # Context/stdres line, e.g.: "  AC>GT  514  145.56  32.40"
  if (!is.na(current_mut) && grepl("^  [ACGT]{2}>[ACGT]{2}\\s", line)) {
    tokens <- strsplit(trimws(line), "\\s+")[[1]]
    records[[length(records) + 1]] <- list(
      section = current_section,
      mutType = current_mut,
      p_value = current_pval,
      context = tokens[1],
      stdres  = suppressWarnings(as.numeric(tokens[4]))
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

# Restrict to target doublet type, drop NA stdres (skip entries)
sub <- df[df$mutType == target_mut & !is.na(df$stdres), ]

if (nrow(sub) == 0) {
  cat(sprintf("No data found for doublet type '%s'.\n", target_mut))
  quit(status = 0)
}

# Rank contexts within each section (rank 1 = highest stdres)
sub$rank <- ave(-sub$stdres, sub$section, FUN = rank)

# Find sections where ALL target contexts rank within top K
valid_sections <- Reduce(intersect, lapply(target_ctxs, function(ctx) {
  sub$section[sub$context == ctx & sub$rank <= K]
}))

cat(sprintf("Filter: %s, contexts = %s (top %d in category)\n\n",
            target_mut, paste(target_ctxs, collapse = ", "), K))

if (length(valid_sections) == 0) {
  cat("No matching sections found.\n")
  quit(status = 0)
}

# Build output: one row per valid section, one stdres column per target context
pval_by_section <- tapply(sub$p_value, sub$section, function(x) x[1])

out <- data.frame(
  section = valid_sections,
  p_value = as.numeric(pval_by_section[valid_sections]),
  stringsAsFactors = FALSE
)
for (ctx in target_ctxs) {
  ctx_sub <- sub[sub$context == ctx, c("section", "stdres")]
  out[[ctx]] <- ctx_sub$stdres[match(out$section, ctx_sub$section)]
}

# Sort by stdres of first specified context (descending)
out <- out[order(out[[target_ctxs[1]]], decreasing = TRUE), ]

sig_col <- ifelse(is.na(out$p_value), "",
           ifelse(out$p_value < 0.01, "**",
           ifelse(out$p_value < 0.05, "*", "")))

# Print header
ctx_header <- paste(sprintf("%8s", target_ctxs), collapse = "  ")
cat(sprintf("%-52s  %12s  %-2s  %s\n", "section", "p_value", "sig", ctx_header))
cat(paste(rep("-", 52 + 16 + 4 + K * 10), collapse = ""), "\n")

for (i in seq_len(nrow(out))) {
  ctx_vals <- paste(sprintf("%8.2f", sapply(target_ctxs, function(ctx) out[[ctx]][i])),
                   collapse = "  ")
  cat(sprintf("%-52s  %12.4e  %-2s  %s\n",
              out$section[i], out$p_value[i], sig_col[i], ctx_vals))
}
cat("\n")
