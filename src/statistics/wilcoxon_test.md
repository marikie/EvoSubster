**2群比較専用**で、しかも **96パターン全部を一括で比較**する R スクリプト案を出せます。
そして結論からいうと：

- **可能です**
- ただし **96個同時検定**になるので、**FDR補正（BH）**は必須です
- 主張は
  **「sampled cnidarian species では sampled fungal species より、その pattern の log2 enrichment が高い」**
  と書くのが科学的に安全です

---

# まず、何を比較するか

各 species ごとに、各 SBS96 pattern について

\[
\text{enrichment}=
\frac{\text{その context の置換率}}
{\text{同じ mutation class の他15 context の置換率}}
\]

を計算します。

例：`ACA>AAA`

\[
\text{rate}_{target} = \frac{\#(ACA>AAA)}{\#(ACA\ original)}
\]

\[
\text{rate}_{bg} = \frac{\sum \#(N[C>A]N)-\#(ACA>AAA)}
{\sum \#(NCN\ original)-\#(ACA\ original)}
\]

\[
\text{enrichment} = \frac{\text{rate}_{target}}{\text{rate}_{bg}}
\]

比較には

\[
\log_2(\text{enrichment})
\]

を使います。

これで species ごとに 96 個の値が出るので、
各 pattern ごとに **group1 vs group2 の Wilcoxon rank-sum test** を行えます。

---

# おすすめの使い方

## もし本当に言いたいのが ACA>AAA / ACG>AAG なら
- まずその **2つを仮説ベースで比較**
- 補助的に **96全部 exploratory に走査**

がいちばんきれいです。

なぜなら、96全部やると多重検定で power が下がるからです。

---

# 2群・96パターン一括比較用 R スクリプト案

## 想定する使い方

```bash
Rscript src/statistics/compare_sbs96_enrichment_2groups.R \
  cnidaria=path/to/cnidaria_tsv_dir \
  fungi=path/to/fungi_tsv_dir
```

- 各ディレクトリには比較したい `.tsv` を入れる
- できれば **全部同じ種類の TSV**
  - 例: 全部 `*_maflinked_ncds.tsv`
- `dinuc.tsv` は自動除外

---

## スクリプト本体

```r
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
```

---

# このスクリプトで何が出るか

## 1. `sbs96_enrichment_per_species.tsv`
各 species × 各 pattern の表です。

重要列：

- `target_label` : 例 `ACA>AAA`
- `mut_type` : 例 `C>A`
- `group`
- `species_id`
- `enrichment`
- `log2_enrichment`

---

## 2. `sbs96_wilcoxon_results.tsv`
各 pattern ごとの群間比較結果です。

重要列：

- `target_label`
- `median_group1`
- `median_group2`
- `median_diff_g1_minus_g2`
- `p_value`
- `p_adj_BH`

### ここで見るべきもの
- `p_adj_BH < 0.05`
- `median_diff_g1_minus_g2 > 0`

なら

> group1 のほうが group2 より、その pattern の enrichment が高い

と解釈できます。

---

# 96パターン全部をまとめて比較することは可能？

## はい、可能です
このスクリプトはそのためのものです。

ただし、**96個まとめて有意差検定**をするときは以下に注意です。

---

## 注意1：多重検定
96個テストするので、**未補正 p値はそのまま使わない**でください。

- おすすめ：**BH (FDR)**
- より厳しくするなら Holm

このスクリプトでは `p_adj_BH` を出しています。

---

## 注意2：生物学的な主張は少し慎重に
この方法は **species を単位**にしているので、かなりまともですが、

- 系統非独立性
- 同じ genus 内 species の近縁性
- 同じ species が複数 TSV で入る問題

は残ります。

だから主張は

> **sampled cnidarian species** では
> **sampled fungal species** より
> ACA>AAA / ACG>AAG の enrichment が高い

と書くのが安全です。

---

## 注意3：同じ species を複数回入れない
これはかなり重要です。

たとえば *Acropora hyacinthus* が複数 trio から入ると、
独立サンプルとして重複カウントになります。

### 推奨
- 1 focal species = 1 TSV にする
- あるいは species ごとに平均してから比較

---

# あなたの主張をこの結果でどう書けるか

## もし ACA>AAA と ACG>AAG が BH補正後も有意なら
こう書けます。

> sampled species を単位とした比較では、ACA>AAA および ACG>AAG の log2 enrichment は cnidaria で fungi より高かった（Wilcoxon rank-sum test, BH-adjusted p < 0.05）。

より自然にするなら：

> sampled cnidarian species では、sampled fungal species と比べて、ACA>AAA および ACG>AAG が C>A class の中で相対的に高頻度であった。

---

## 有意でなければ
> 複数の cnidarian species で ACA>AAA および ACG>AAG の上昇は繰り返し観察されたが、species-level の群 間比較では fungi との差は明確ではなかった。したがって現時点では、これらを cnidaria-specific と断定するよりも、cnidaria で繰り返し見られる候補的パターンとして扱うのが適切である。

---

# 実務上のおすすめ

あなたの目的なら、順番としては：

1. **まず ACA>AAA / ACG>AAG の2つだけ比較**
2. 次に **96全部を exploratory に走査**
3. 補助図として top hit を boxplot

がいちばんきれいです。

---

# 3群以上は？
今回のスクリプトは 2群専用です。
3群以上は前に書いた通り、

- **Kruskal–Wallis**
- その後 **pairwise Wilcoxon**

が簡単でよいです。

---

必要なら次に、

1. **このスクリプトを “target 2パターンだけ版” に簡略化**
2. **TSVの実際の列名に合わせて repo 用に微修正版を書く**
3. **結果を Figure 用に並べる ggplot 版を書く**

のどれかをすぐできます。

