library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(showtext)
library(sysfonts)
library(rlang)

generate_plots <- function(tsv_path) {
  data <- read.csv(tsv_path, sep = "\t", header = TRUE)

  path_without_extension <- tools::file_path_sans_ext(tsv_path)
  # graph_path_exp <- paste(path_without_extension, "_logRatio_exp.pdf", sep = "")
  graph_path_mean <- paste(path_without_extension, "_logRatio.pdf", sep = "")
  data <- add_logRatio(data)
  # create_pdf(graph_path_exp, data, value_col = "logRatio_exp")
  create_pdf(graph_path_mean, data, value_col = "logRatio_mean")
}

create_pdf <- function(graph_path, data, value_col) {
  # Prepare ordering by substitution type groups and within-group order
  group_order <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  # Fonts for labels
  sysfonts::font_add_google("Courier Prime", "mn", 700)
  sysfonts::font_add_google("Roboto", "os")
  showtext::showtext_auto()
  data <- data %>%
    mutate(trans = factor(substr(mutType, 3, 5), levels = group_order)) %>%
    arrange(trans, oriType)

  # Insert small gaps between groups on the x-axis
  gap_width <- 1
  data <- data %>%
    group_by(trans) %>%
    mutate(pos_in_group = row_number()) %>%
    ungroup()

  group_info <- data %>%
    group_by(trans) %>%
    summarise(group_size = max(pos_in_group), .groups = "drop") %>%
    arrange(trans) %>%
    mutate(offset = cumsum(dplyr::lag(group_size + gap_width, default = 0)))

  data <- data %>%
    left_join(group_info, by = "trans") %>%
    mutate(pos = pos_in_group + offset)

  value_sym <- rlang::sym(value_col)
  data[[value_col]][!is.finite(data[[value_col]])] <- NA_real_
  finite_values <- data[[value_col]][is.finite(data[[value_col]])]
  if (length(finite_values) == 0) {
    warning(sprintf("No finite values detected for %s; defaulting to zeros.", value_col))
    finite_values <- 0
  }

  # X labels: oriType in the arranged order
  x_breaks <- data$pos
  x_labels <- data$oriType

  # Header rectangles like sbmut.R's add_colored_rectangles
  color_array <- RColorBrewer::brewer.pal(6, "Dark2")
  fill_map <- setNames(color_array, group_order)

  yr <- range(finite_values, na.rm = TRUE)
  max_abs <- max(abs(yr))
  if (!is.finite(max_abs) || max_abs == 0) {
    max_abs <- 1
  }
  # Create symmetric range around 0
  ylim_high <- max_abs * 1.4 # Space above for rectangles and text
  ylim_low <- -max_abs * 1.4 # Equal space below for labels and grid

  # Calculate max integer for grid breaks
  max_int <- max(1, ceiling(max_abs))

  # Position rectangles and labels relative to symmetric limits
  yrect_low <- ylim_high - 0.10 * (ylim_high - ylim_low)
  yrect_high <- ylim_high - 0.06 * (ylim_high - ylim_low)
  y_text <- ylim_high - 0.02 * (ylim_high - ylim_low)
  y_label3 <- ylim_low + 0.02 * (ylim_high - ylim_low)
  y_label2 <- ylim_low - 0.01 * (ylim_high - ylim_low)
  y_label1 <- ylim_low - 0.04 * (ylim_high - ylim_low)

  gdf <- data %>%
    group_by(trans) %>%
    summarise(
      xmin = min(pos) - 0.5,
      xmax = max(pos) + 0.5,
      xmid = (min(pos) + max(pos)) / 2,
      .groups = "drop"
    ) %>%
    mutate(trans = factor(trans, levels = group_order)) %>%
    arrange(trans) %>%
    mutate(label = c("C > A", "C > G", "C > T", "T > A", "T > C", "T > G"))

  labels_df <- data %>%
    distinct(pos, oriType, trans) %>%
    mutate(
      c1 = substr(oriType, 1, 1),
      c2 = substr(oriType, 2, 2),
      c3 = substr(oriType, 3, 3)
    )
  if (value_col == "logRatio_exp") {
    y_label <- "Log2{(#sbst/#ori) / (expected #sbst/#ori)}"
  }
  if (value_col == "logRatio_mean") {
    y_label <- "Log2{(#sbst/#ori) / (mean of #sbst/#ori)}"
  }
  p <- ggplot(data, aes(x = pos, y = !!value_sym)) +
    geom_rect(
      data = gdf,
      aes(xmin = xmin, xmax = xmax, ymin = yrect_low, ymax = yrect_high, fill = trans),
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    scale_fill_manual(values = fill_map, guide = "none") +
    # Make the header labels bold by setting fontface = "bold"
    annotate("text", x = gdf$xmid, y = y_text, label = gdf$label, size = 10, family = "os", fontface = "bold") +
    geom_segment(
      data = data,
      aes(x = pos, xend = pos, y = 0, yend = !!value_sym, color = trans),
      inherit.aes = FALSE,
      linewidth = 0.8, alpha = 0.9, show.legend = FALSE, na.rm = TRUE
    ) +
    scale_color_manual(values = fill_map, guide = "none") +
    geom_point(size = 4, aes(color = trans), show.legend = FALSE, na.rm = TRUE) +
    # Custom x labels: outer letters in black, middle letter in group color
    geom_text(
      data = labels_df, aes(x = pos, y = y_label3, label = c3),
      inherit.aes = FALSE, angle = 90, size = 8, family = "mn", color = "black"
    ) +
    geom_text(
      data = labels_df, aes(x = pos, y = y_label1, label = c1),
      inherit.aes = FALSE, angle = 90, size = 8, family = "mn", color = "black"
    ) +
    geom_text(
      data = subset(labels_df, trans == "C>A"), aes(x = pos, y = y_label2, label = c2),
      inherit.aes = FALSE, angle = 90, size = 8, family = "mn", color = fill_map["C>A"]
    ) +
    geom_text(
      data = subset(labels_df, trans == "C>G"), aes(x = pos, y = y_label2, label = c2),
      inherit.aes = FALSE, angle = 90, size = 8, family = "mn", color = fill_map["C>G"]
    ) +
    geom_text(
      data = subset(labels_df, trans == "C>T"), aes(x = pos, y = y_label2, label = c2),
      inherit.aes = FALSE, angle = 90, size = 8, family = "mn", color = fill_map["C>T"]
    ) +
    geom_text(
      data = subset(labels_df, trans == "T>A"), aes(x = pos, y = y_label2, label = c2),
      inherit.aes = FALSE, angle = 90, size = 8, family = "mn", color = fill_map["T>A"]
    ) +
    geom_text(
      data = subset(labels_df, trans == "T>C"), aes(x = pos, y = y_label2, label = c2),
      inherit.aes = FALSE, angle = 90, size = 8, family = "mn", color = fill_map["T>C"]
    ) +
    geom_text(
      data = subset(labels_df, trans == "T>G"), aes(x = pos, y = y_label2, label = c2),
      inherit.aes = FALSE, angle = 90, size = 8, family = "mn", color = fill_map["T>G"]
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "#666666", linewidth = 0.6) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0, 0), minor_breaks = NULL) +
    scale_y_continuous(breaks = seq(-max_int, max_int, by = 1), minor_breaks = NULL) +
    labs(x = "Original Trinucleotides", y = y_label) +
    theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 26, margin = margin(t = 10), family = "os"),
      axis.title.y = element_text(size = 26, margin = margin(r = 10), family = "os"),
      axis.text.y = element_text(size = 22, margin = margin(r = 8)),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(color = "#cecece", linewidth = 0.2),
      panel.grid.major.y = element_line(color = "#cecece", linewidth = 0.5),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(t = 12, r = 12, b = 12, l = 20)
    ) +
    coord_cartesian(ylim = c(ylim_low, ylim_high), clip = "off")

  ggsave(filename = graph_path, plot = p, device = "pdf", width = 30, height = 8)
}

add_logRatio <- function(data) {
  # Add oriType column to data
  data <- data %>%
    mutate(oriType = paste0(substr(mutType, 1, 1), substr(mutType, 3, 3), substr(mutType, nchar(mutType), nchar(mutType))))

  # Add transition column to data
  # Add observed figure of mutNum/totalRootNum
  obs_mut_over_ori <- ifelse(is.finite(data$totalRootNum) & data$totalRootNum > 0,
    data$mutNum / data$totalRootNum,
    NA_real_
  )

  # all_sbst_sum <- data %>%
  #   select(mutNum) %>%
  #   sum()
  # all_ori_sum <- data %>%
  #   select(oriType, totalRootNum) %>%
  #   unique() %>%
  #   pull(totalRootNum) %>%
  #   sum()
  # expected_sbst_over_ori <- (all_sbst_sum / all_ori_sum) / 3
  positive_rates <- obs_mut_over_ori[is.finite(obs_mut_over_ori) & obs_mut_over_ori > 0]
  mean_obs_sbst_over_ori <- if (length(positive_rates) > 0) mean(positive_rates) else NA_real_
  if (!is.finite(mean_obs_sbst_over_ori) || mean_obs_sbst_over_ori <= 0) {
    warning("Mean observed substitution rate is non-positive; log ratios will be NA.")
    mean_obs_sbst_over_ori <- NA_real_
  }

  if (!is.na(mean_obs_sbst_over_ori) && mean_obs_sbst_over_ori > 0) {
    ratio_base <- obs_mut_over_ori / mean_obs_sbst_over_ori
    ratio_base[!is.finite(ratio_base) | ratio_base <= 0] <- NA_real_
    log_ratio_mean <- log2(ratio_base)
  } else {
    log_ratio_mean <- rep(NA_real_, length(obs_mut_over_ori))
  }

  data <- data %>%
    mutate(
      obs_mut_over_ori = obs_mut_over_ori,
      # logRatio_exp = log2(obs_mut_over_ori / expected_sbst_over_ori),
      logRatio_mean = log_ratio_mean
    )
  # print(paste("logRatio_exp: ", data$logRatio_exp))
  # print(paste("logRatio_mean: ", data$logRatio_mean))
  # print(paste("expected_sbst_over_ori: ", expected_sbst_over_ori))
  # print(paste("mean_obs_sbst_over_ori: ", mean_obs_sbst_over_ori))
  return(data)
}

# Access the arguments
args <- commandArgs(trailingOnly = TRUE)

# Access the arguments
tsv_path <- args[1] # File path for the input data, .tsv file
generate_plots(tsv_path)
