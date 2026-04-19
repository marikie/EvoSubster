library(stringr)
library(RColorBrewer)
library(showtext)
library(sysfonts)
library(dplyr)

# MODIFIED THE CODE FROM https://github.com/kartong88/Plot-Mutation-Landscape

generate_plots <- function(tsv_path, out_dir = dirname(tsv_path)) {
  data <- read_and_transform_data(tsv_path)
  trinuc.lab <- data %>% pull(oriType)

  # out_dir is the singlenuc parent; split outputs into ratio/ and count/
  ratio_dir <- file.path(out_dir, "ratio")
  count_dir <- file.path(out_dir, "count")
  dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(count_dir, recursive = TRUE, showWarnings = FALSE)
  base <- tools::file_path_sans_ext(basename(tsv_path))
  normgraph_path <- file.path(ratio_dir, paste0(base, "_norm.pdf"))
  sbstgraph_path <- file.path(count_dir, paste0(base, "_sbst.pdf"))
  origraph_path  <- file.path(count_dir, paste0(base, "_ori.pdf"))

  print(paste("Creating norm graph...", normgraph_path))
  create_pdf_plot(normgraph_path, data, trinuc.lab, "norm")
  print(paste("Creating sbst graph...", sbstgraph_path))
  create_pdf_plot(sbstgraph_path, data, trinuc.lab, "sbst")
  print(paste("Creating ori graph...", origraph_path))
  create_pdf_plot_ori(origraph_path, data, trinuc.lab)

  process_norm_graph(data)
}

read_and_transform_data <- function(file_path) {
  data <- read.csv(file_path, sep = "\t", header = TRUE) %>%
    as_tibble()

  # Add oriType column using string manipulation on mutType
  data <- data %>%
    mutate(
      oriType = paste0(
        str_sub(mutType, 1, 1),
        str_sub(mutType, 3, 3),
        str_sub(mutType, 7, 7)
      )
    )

  return(data)
}

extract_values <- function(data, graph_type) {
  if (graph_type == "norm") {
    return(data %>%
      mutate(value = mutNum / totalRootNum * 100) %>%
      pull(value))
  } else if (graph_type == "sbst") {
    return(data %>% pull(mutNum))
  }
}

compute_y_axis_max <- function(conv_data_norm) {
  raw_max <- max(conv_data_norm, na.rm = TRUE)
  if (ceiling(raw_max) == 1) {
    if (raw_max < 0.1) {
      return(round(raw_max, 2))
    } else {
      return(round(raw_max, 1))
    }
  } else {
    return(ceiling(raw_max))
  }
}

create_pdf_plot <- function(filename, data, trinuc.lab, graph_type) {
  conv.data.norm <- extract_values(data, graph_type)
  pct_yaxs_max <- compute_y_axis_max(conv.data.norm)
  color_array <- brewer.pal(6, "Dark2")
  pdf(filename, width = 30, height = 8)
  font_add_google("Courier Prime", "mn", 700)
  font_add_google("Roboto", "os")
  showtext_auto()
  par(family = "mn", mar = c(7.5, 6, 4, 1), cex.axis = 2)
  bar_positions <- barplot(
    axes = FALSE,
    family = "os",
    as.numeric(conv.data.norm),
    col = rep(color_array, each = 16),
    cex.names = 0.7,
    las = 3,
    ylim = c(0, pct_yaxs_max + 0.15 * pct_yaxs_max),
    space = 0,
  )
  # manually add y axis
  axis(side = 2, line = -3.5)
  add_labels(trinuc.lab, bar_positions)
  add_axis_names(graph_type)
  add_colored_rectangles(bar_positions, pct_yaxs_max, color_array)
  dev.off()
}

create_pdf_plot_ori <- function(filename, data, trinuc.lab) {
  # Get unique oriType and totalRootNum combinations
  unique_data <- data %>%
    select(oriType, totalRootNum) %>%
    unique()

  unique_labels <- unique_data %>% pull(oriType)
  aggregated_counts <- unique_data %>% pull(totalRootNum)
  pct_yaxs_max <- compute_y_axis_max(aggregated_counts)
  upper_limit <- if (pct_yaxs_max == 0) 1 else pct_yaxs_max + 0.15 * pct_yaxs_max
  pdf(filename, width = 30 * (32 / 96), height = 8)
  font_add_google("Courier Prime", "mn", 700)
  font_add_google("Roboto", "os")
  showtext_auto()
  par(family = "mn", mar = c(7.5, 6, 4, 1), cex.axis = 2)
  bar_positions <- barplot(
    axes = FALSE,
    family = "os",
    as.numeric(aggregated_counts),
    cex.names = 1.2,
    col = "#b4b4b4",
    border = "black",
    names.arg = rep("", length(unique_labels)),
    ylim = c(0, upper_limit),
    space = 0,
  )
  axis(side = 2, line = -0.5)
  add_unique_labels(unique_labels, bar_positions)
  add_axis_names("ori")
  dev.off()
}

add_labels <- function(trinuc.lab, bar_positions) {
  for (i in seq_along(trinuc.lab)) {
    label <- trinuc.lab[i]
    x_pos <- bar_positions[i]
    base_y <- par("usr")[3] - 0.03 * (par("usr")[4] - par("usr")[3])
    if (i %% 16 != 0) {
      coln <- i %/% 16 + 1
    } else {
      coln <- i %/% 16
    }
    text(family = "mn", x_pos, base_y, substr(label, 3, 3), xpd = TRUE, cex = 2, srt = 90)
    text(family = "mn", x_pos, base_y - strwidth(label) * 0.03 * (par("usr")[4] - par("usr")[3]), substr(label, 2, 2), xpd = TRUE, cex = 2, col = brewer.pal(6, "Dark2")[coln], srt = 90)
    text(family = "mn", x_pos, base_y - strwidth(label) * 0.06 * (par("usr")[4] - par("usr")[3]), substr(label, 1, 1), xpd = TRUE, cex = 2, srt = 90)
  }
}

add_unique_labels <- function(labels, bar_positions) {
  if (length(labels) == 0) {
    return()
  }
  y_range <- par("usr")[4] - par("usr")[3]
  base_y <- par("usr")[3] - 0.07 * y_range
  text(family = "mn", x = bar_positions, y = base_y, labels = labels, xpd = TRUE, cex = 2, srt = 90)
}

add_axis_names <- function(graph_type) {
  if (graph_type == "norm") {
    mtext(family = "os", "Original Trinucleotides", side = 1, line = 5, cex = 2.5)
    mtext(family = "os", "#Sbst/#OriTrinucs (%)", side = 2, line = 0.5, cex = 2.5)
  } else if (graph_type == "sbst") {
    mtext(family = "os", "Original Trinucleotides", side = 1, line = 5, cex = 2.5)
    mtext(family = "os", "#Substitutions", side = 2, line = 0.5, cex = 2.5)
  } else if (graph_type == "ori") {
    mtext(family = "os", "Original Trinucleotide Patterns", side = 1, line = 5, cex = 2.5)
    mtext(family = "os", "#Original Trinucleotides", side = 2, line = 2.5, cex = 2.5)
  }
}

add_colored_rectangles <- function(bar_positions, pct_yaxs_max, color_array) {
  total_size_per_group <- mean(diff(bar_positions)) * 16
  for (i in 1:6) {
    left <- (i - 1) * total_size_per_group + 0.2
    right <- i * total_size_per_group - 0.2
    label_mid <- total_size_per_group / 2 + (i - 1) * total_size_per_group
    rect(left, pct_yaxs_max + 0.05 * pct_yaxs_max, right, pct_yaxs_max + 0.08 * pct_yaxs_max, col = color_array[i], border = NA)
    text(family = "mn", x = label_mid, y = pct_yaxs_max + 0.12 * pct_yaxs_max, labels = c("C > A", "C > G", "C > T", "T > A", "T > C", "T > G")[i], cex = 2.5)
  }
}

process_norm_graph <- function(data) {
  # Calculate mutation percentage and create display names
  data <- data %>%
    mutate(
      MutationPercentage = mutNum / totalRootNum * 100,
      display_name = case_when(
        nchar(mutType) >= 7 ~ paste0(
          str_sub(mutType, 1, 1), str_sub(mutType, 3, 3), str_sub(mutType, 7, 7),
          " → ",
          str_sub(mutType, 1, 1), str_sub(mutType, 5, 5), str_sub(mutType, 7, 7)
        ),
        TRUE ~ mutType
      )
    )

  # Sort data by mutation percentage (descending)
  data_sorted <- data %>%
    arrange(desc(MutationPercentage))

  print("Data sorted by Mutation Percentage:")
  print(data_sorted)

  cat("\n\nTop 10:\n")
  data_sorted %>%
    slice_head(n = 10) %>%
    print()

  cat("\n\nWorst 10:\n")
  data_sorted %>%
    slice_tail(n = 10) %>%
    print()

  # Reverse sort
  data_rev_sorted <- data %>%
    arrange(MutationPercentage)

  cat("\n\nWorst 10 (rev):\n")
  data_rev_sorted %>%
    slice_head(n = 10) %>%
    print()
}


args <- commandArgs(trailingOnly = TRUE)

# Access the arguments
tsv_path <- args[1] # File path for the input data, .tsv file
out_dir <- if (length(args) >= 2) args[2] else dirname(tsv_path)
generate_plots(tsv_path, out_dir)
