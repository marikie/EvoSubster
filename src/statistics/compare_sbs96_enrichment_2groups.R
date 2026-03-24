#!/usr/bin/env Rscript

# Compare SBS96 enrichment between two groups of EvoSubster TSV files
#
# Usage:
#   Rscript compare_sbs96_enrichment_2groups.R group1=DIR1 group2=DIR2
#
# Example:
#   Rscript compare_sbs96_enrichment_2groups.R \
#     cnidaria=./cnidaria_tsv \
#     fungi=./fungi_tsv
#
# Outputs:
#   sbs96_enrichment_per_species.tsv
#   sbs96_wilcoxon_results.tsv
#   sbs96_boxplot_top12.pdf
#
# Notes:
# - Recommended: one focal species = one TSV
# - Recommended: use the same TSV type across both groups
# - This script excludes dinuc TSVs automatically

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
"Usage:
  Rscript compare_sbs96_enrichment_2groups.R group1=DIR1 group2=DIR2

Example:
  Rscript compare_sbs96_enrichment_2groups.R \
    cnidaria=./cnidaria_tsv \
    fungi=./fungi_tsv
", file = stderr())
  quit(status = 1)
}

if (length(args) != 2) usage()
if (!all(grepl("=", args, fixed = TRUE))) usage()

specs <- strsplit(args, "=", fixed = TRUE)
group_names <- vapply(specs, `[`, character(1), 1)
group_dirs  <- vapply(specs, `[`, character(1), 2)

if (length(unique(group_names)) != 2) {
  stop("Please provide exactly two distinct group names.")
}
if (any(!dir.exists(group_dirs))) {
  stop("Directory not found: ", paste(group_dirs[!dir.exists(group_dirs)], collapse = ", "))
}

# ----------------------------
# helper functions
# ----------------------------

normalize_label <- function(x) {
  x <- trimws(as.character(x))
  out <- x

  # Convert bracket notation A[C>A]A -> ACA>AAA
  is_bracket <- grepl("^[ACGT]\\[[CT]>[ACGT]\\][ACGT]$", x)
  if (any(is_bracket)) {
    y <- x[is_bracket]
    left  <- substr(y, 1, 1)
    fromb <- substr(y, 3, 3)
    tob   <- substr(y, 5, 5)
    right <- substr(y, 7, 7)
    out[is_bracket] <- paste0(left, fromb, right, ">", left, tob, right)
  }
  out
}

detect_label_col <- function(df) {
  pats <- c("^[ACGT]{3}>[ACGT]{3}$", "^[ACGT]\\[[CT]>[ACGT]\\][ACGT]$")
  best_col <- NA_character_
  best_score <- -1

  for (nm in names(df)) {
    vals <- as.character(df[[nm]])
    score <- mean(Reduce(`|`, lapply(pats, grepl, x = vals)), na.rm = TRUE)
    if (!is.na(score) && score > best_score) {
      best_score <- score
      best_col <- nm
    }
  }

  if (is.na(best_col) || best_score < 0.5) return(NA_character_)
  best_col
}

detect_count_cols <- function(df, label_col) {
  nm_low <- tolower(names(df))
  is_num <- vapply(df, is.numeric, logical(1))

  obs_idx <- grep("sub|mut|count", nm_low)
  bg_idx  <- grep("orig|ori|original|opport|background", nm_low)

  obs_idx <- obs_idx[names(df)[obs_idx] != label_col]
  bg_idx  <- bg_idx[names(df)[bg_idx] != label_col]

  obs_idx <- obs_idx[is_num[obs_idx]]
  bg_idx  <- bg_idx[is_num[bg_idx]]

  obs_col <- if (length(obs_idx) > 0) names(df)[obs_idx[1]] else NA_character_
  bg_col  <- if (length(bg_idx) > 0) names(df)[bg_idx[1]] else NA_character_

  # fallback: first two numeric cols except label
  if (is.na(obs_col) || is.na(bg_col) || obs_col == bg_col) {
    num_cols <- names(df)[is_num]
    num_cols <- setdiff(num_cols, label_col)
    if (length(num_cols) < 2) {
      stop("Could not detect observed/original count columns.")
    }
    obs_col <- num_cols[1]
    bg_col  <- num_cols[2]
  }

  c(obs_col = obs_col, bg_col = bg_col)
}

read_evosubster_tsv <- function(path) {
  df <- tryCatch(
    read.delim(path, header = TRUE, sep = "\t", check.names = FALSE,
               stringsAsFactors = FALSE, comment.char = "", quote = ""),
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) < 2) {
    stop("Failed to read TSV: ", path)
  }

  label_col <- detect_label_col(df)

  # fallback: header = FALSE
  if (is.na(label_col)) {
    df2 <- tryCatch(
      read.delim(path, header = FALSE, sep = "\t", check.names = FALSE,
                 stringsAsFactors = FALSE, comment.char = "", quote = ""),
      error = function(e) NULL
    )
    if (!is.null(df2)) {
      names(df2) <- paste0("V", seq_len(ncol(df2)))
      label_col2 <- detect_label_col(df2)
      if (!is.na(label_col2)) {
        df <- df2
        label_col <- label_col2
      }
    }
  }

  if (is.na(label_col)) {
    stop("Could not detect context label column in: ", path)
  }

  count_cols <- detect_count_cols(df, label_col)

  out <- data.frame(
    label = normalize_label(df[[label_col]]),
    obs_count = suppressWarnings(as.numeric(df[[count_cols["obs_col"]]])),
    orig_count = suppressWarnings(as.numeric(df[[count_cols["bg_col"]]])),
    stringsAsFactors = FALSE
  )

  out <- out[grepl("^[ACGT]{3}>[ACGT]{3}$", out$label), , drop = FALSE]
  if (nrow(out) == 0) {
    stop("No valid SBS rows found in: ", path)
  }

  out$anc <- substr(out$label, 1, 3)
  out$der <- substr(out$label, 5, 7)
  out$mut_type <- paste0(substr(out$anc, 2, 2), ">", substr(out$der, 2, 2))

  out
}

extract_species_id <- function(path) {
  sub("\\.tsv$", "", basename(path), ignore.case = TRUE)
}

compute_all96_enrichment <- function(df) {
  # Haldane-style smoothing to avoid Inf / NA for very sparse patterns
  out_list <- vector("list", nrow(df))

  for (i in seq_len(nrow(df))) {
    target_label <- df$label[i]
    target_mut_type <- df$mut_type[i]

    class_rows <- df[df$mut_type == target_mut_type, , drop = FALSE]

    obs_target <- df$obs_count[i]
    opp_target <- df$orig_count[i]

    obs_bg <- sum(class_rows$obs_count, na.rm = TRUE) - obs_target
    opp_bg <- sum(class_rows$orig_count, na.rm = TRUE) - opp_target

    # Smoothed rates
    rate_target <- (obs_target + 0.5) / (opp_target + 1)
    rate_bg     <- (obs_bg + 0.5) / (opp_bg + 1)

    enrichment <- rate_target / rate_bg
    log2_enrichment <- log2(enrichment)

    out_list[[i]] <- data.frame(
      target_label = target_label,
      mut_type = target_mut_type,
      anc = df$anc[i],
      der = df$der[i],
      target_obs = obs_target,
      target_opp = opp_target,
      bg_obs_other15 = obs_bg,
      bg_opp_other15 = opp_bg,
      target_rate = rate_target,
      bg_rate_other15 = rate_bg,
      enrichment = enrichment,
      log2_enrichment = log2_enrichment,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out_list)
}

# ----------------------------
# collect all per-species enrichment
# ----------------------------

all_res <- list()
k <- 1

for (g in seq_along(group_names)) {
  grp <- group_names[g]
  dirp <- group_dirs[g]

  files <- list.files(dirp, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("dinuc", basename(files), ignore.case = TRUE)]

  if (length(files) == 0) {
    warning("No TSV files found in ", dirp)
    next
  }

  for (f in files) {
    dat <- tryCatch(read_evosubster_tsv(f), error = function(e) {
      warning("Skipping ", f, " : ", conditionMessage(e))
      NULL
    })
    if (is.null(dat)) next

    enr <- compute_all96_enrichment(dat)
    enr$group <- grp
    enr$file <- f
    enr$species_id <- extract_species_id(f)

    all_res[[k]] <- enr
    k <- k + 1
  }
}

if (length(all_res) == 0) {
  stop("No valid data collected.")
}

res <- do.call(rbind, all_res)

write.table(
  res,
  file = "sbs96_enrichment_per_species.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ----------------------------
# Wilcoxon test for each SBS96 pattern
# ----------------------------

g1 <- group_names[1]
g2 <- group_names[2]

contexts <- sort(unique(res$target_label))
test_rows <- list()
k <- 1

for (ctx in contexts) {
  sub <- res[res$target_label == ctx & is.finite(res$log2_enrichment), , drop = FALSE]

  x1 <- sub$log2_enrichment[sub$group == g1]
  x2 <- sub$log2_enrichment[sub$group == g2]

  n1 <- length(x1)
  n2 <- length(x2)

  # Require at least 2 species per group for a meaningful comparison
  if (n1 < 2 || n2 < 2) {
    test_rows[[k]] <- data.frame(
      target_label = ctx,
      mut_type = unique(sub$mut_type)[1],
      n_group1 = n1,
      n_group2 = n2,
      median_group1 = if (n1 > 0) median(x1) else NA_real_,
      median_group2 = if (n2 > 0) median(x2) else NA_real_,
      median_diff_g1_minus_g2 = if (n1 > 0 && n2 > 0) median(x1) - median(x2) else NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    )
    k <- k + 1
    next
  }

  wt <- wilcox.test(x1, x2, exact = FALSE)

  test_rows[[k]] <- data.frame(
    target_label = ctx,
    mut_type = unique(sub$mut_type)[1],
    n_group1 = n1,
    n_group2 = n2,
    median_group1 = median(x1),
    median_group2 = median(x2),
    median_diff_g1_minus_g2 = median(x1) - median(x2),
    statistic = unname(wt$statistic),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
  k <- k + 1
}

tests <- do.call(rbind, test_rows)
tests$p_adj_BH <- p.adjust(tests$p_value, method = "BH")
tests <- tests[order(tests$p_adj_BH, tests$p_value), ]

write.table(
  tests,
  file = "sbs96_wilcoxon_results.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ----------------------------
# simple plot: top 12 patterns by BH-adjusted p
# ----------------------------

plot_table <- tests[is.finite(tests$p_adj_BH), , drop = FALSE]
plot_table <- plot_table[order(plot_table$p_adj_BH, plot_table$p_value), , drop = FALSE]
topN <- min(12, nrow(plot_table))

if (topN > 0) {
  top_contexts <- plot_table$target_label[seq_len(topN)]
  pdf("sbs96_boxplot_top12.pdf", width = 10, height = 3 * topN)
  par(mfrow = c(topN, 1), mar = c(5, 4, 2.5, 1))

  for (ctx in top_contexts) {
    sub <- res[res$target_label == ctx & is.finite(res$log2_enrichment), , drop = FALSE]
    boxplot(log2_enrichment ~ group, data = sub,
            main = paste0(ctx, " (", unique(sub$mut_type), ")"),
            ylab = "log2 enrichment",
            xlab = "",
            outline = FALSE,
            col = "grey90",
            border = "grey40")
    stripchart(log2_enrichment ~ group, data = sub,
               method = "jitter", vertical = TRUE,
               pch = 16, col = "#3366AA88", add = TRUE)
    abline(h = 0, lty = 2, col = "red")
  }
  dev.off()
}

cat("Done.\n")
cat("Output files:\n")
cat("  sbs96_enrichment_per_species.tsv\n")
cat("  sbs96_wilcoxon_results.tsv\n")
cat("  sbs96_boxplot_top12.pdf\n")