# Helpers loaded by trio_selection.R into its calling environment.

accession_tip_labels <- function(tree, leaves) {
  rows <- match(tree$tip.label, leaves$tip)
  if (anyNA(rows)) {
    stop("No matching row in leaves for tree tip(s): ",
         paste(tree$tip.label[is.na(rows)], collapse = ", "), call. = FALSE)
  }
  paste0(tree$tip.label, "_", leaves$accession[rows])
}

tree_pdf_size <- function(ntip) {
  c(width = 11, height = max(7, 0.18 * ntip + 2))
}

# grDevices::pdf() leaves Helvetica non-embedded, so some viewers substitute a
# poorly rendered fallback that looks blurry when zoomed. cairo_pdf() embeds a
# font subset for consistent, crisp text; fontconfig maps absent Helvetica to a
# metric-compatible sans such as Nimbus Sans.
open_tree_pdf <- function(path, size) {
  if (!capabilities("cairo")) {
    stop("R was built without cairo; cannot write the tree PDF with embedded fonts: ", path, call. = FALSE)
  }
  grDevices::cairo_pdf(path, width = size[["width"]], height = size[["height"]], family = "Helvetica")
}

write_stage0_tree <- function(tree, leaves, out_dir) {
  labeled <- tree
  labeled$tip.label <- accession_tip_labels(tree, leaves)
  files <- c(nwk = file.path(out_dir, "pruned_tree.nwk"),
             pdf = file.path(out_dir, "pruned_tree.pdf"))
  ape::write.tree(labeled, file = files[["nwk"]])

  ntip <- ape::Ntip(labeled)
  size <- tree_pdf_size(ntip)
  open_tree_pdf(files[["pdf"]], size)
  on.exit(grDevices::dev.off(), add = TRUE)
  ape::plot.phylo(labeled, cex = 0.7, label.offset = 0.01, underscore = TRUE)
  title(sprintf("Stage 0: %d species, one assembly each", ntip))
  invisible(files)
}

trio_tip_roles <- function(selected) {
  if (is.null(nrow(selected)) || nrow(selected) == 0L) {
    return(data.frame(tip = character(), trio = integer(), role = character(),
                      stringsAsFactors = FALSE))
  }
  data.frame(
    tip = as.vector(t(as.matrix(selected[, c("out_tip", "in1_tip", "in2_tip")]))),
    trio = rep(seq_len(nrow(selected)), each = 3L),
    role = rep(c("out", "in1", "in2"), times = nrow(selected)),
    stringsAsFactors = FALSE
  )
}

stage1_tree_paths <- function(out_dir) {
  c(nwk = file.path(out_dir, "selected_trios.nwk"),
    pdf = file.path(out_dir, "selected_trios.pdf"))
}

remove_stage1_tree_outputs <- function(out_dir) {
  files <- stage1_tree_paths(out_dir)
  unlink(files)
  invisible(files)
}

selected_tree_pdf_layout <- function(tree, tip_idx, tip_text, n_trio, execution_info = NULL) {
  ntip <- ape::Ntip(tree)
  size <- tree_pdf_size(ntip)
  # Measure with the same embedded font and text sizes as the final PDF.
  measure_pdf <- tempfile("trio-tree-measure-", fileext = ".pdf")
  open_tree_pdf(measure_pdf, size)
  on.exit({
    grDevices::dev.off()
    unlink(measure_pdf)
  }, add = TRUE)

  marker_steps <- 0.22 + strwidth(tip_text, units = "inches", cex = 0.6)
  tip_widths <- numeric(ntip)
  for (i in seq_along(tip_idx)) {
    tip_widths[tip_idx[i]] <- tip_widths[tip_idx[i]] + marker_steps[i]
  }
  annotation_width <- max(tip_widths) + 0.1
  label_width <- max(strwidth(tree$tip.label, units = "inches", cex = 0.7, font = 3))

  # Grow the page by the extra annotation width, preserving the branch space
  # from the original T<number> layout (and at least one inch for long labels).
  old_marker_step <- 0.22 + max(c(0, strwidth(paste0("T", seq_len(n_trio)),
                                            units = "inches", cex = 0.6)))
  old_annotation_width <- max(tabulate(tip_idx, nbins = ntip)) * old_marker_step + 0.1
  size[["width"]] <- size[["width"]] + max(
    0, annotation_width - old_annotation_width,
    annotation_width + label_width + 1 - par("pin")[[1]]
  )
  # Reserve physical space for every metadata line, without shrinking the tree.
  execution_lines <- wrap_execution_lines(execution_info_lines(execution_info), size[["width"]] - 1.2)
  line_height <- 0.16
  header_extra <- 0.8 + length(execution_lines) * line_height
  size[["height"]] <- size[["height"]] + header_extra
  list(size = size, marker_steps = marker_steps, header_extra = header_extra,
       annotation_width = annotation_width, label_width = label_width,
       execution_lines = execution_lines, line_height = line_height)
}

draw_selected_tree_header <- function(n_trio, layout) {
  old_xpd <- par(xpd = NA)
  on.exit(par(old_xpd), add = TRUE)
  usr <- par("usr")
  x_per_inch <- diff(grconvertX(c(0, 1), from = "inches", to = "user"))
  y_per_inch <- diff(grconvertY(c(0, 1), from = "inches", to = "user"))
  center <- mean(usr[1:2])
  metadata_height <- length(layout$execution_lines) * layout$line_height
  text(center, usr[4] + (1.35 + metadata_height) * y_per_inch,
       sprintf("Stage 1: %d %s selected", n_trio, if (n_trio == 1L) "trio" else "trios"),
       cex = 1.2, font = 2)
  text_x <- grconvertX(0.6, from = "inches", to = "user")
  for (i in seq_along(layout$execution_lines)) {
    y <- usr[4] + (1.05 + (length(layout$execution_lines) - i) * layout$line_height) * y_per_inch
    text(text_x, y, layout$execution_lines[i], adj = c(0, 0.5), cex = 0.72)
  }

  key_labels <- c("Outgroup (T#)",
                  "Ingroup: T#.1 = in1 (genome2); T#.2 = in2 (genome3)")
  key_widths <- strwidth(key_labels, units = "inches", cex = 0.8) + 0.45
  left <- center - sum(key_widths) * x_per_inch / 2
  key_y <- usr[4] + 0.65 * y_per_inch
  for (i in seq_along(key_labels)) {
    marker_x <- left + 0.1 * x_per_inch
    points(marker_x, key_y, pch = 8, cex = 1.6, lwd = 2)
    if (i == 1L) points(marker_x, key_y, pch = 1, cex = 1.6, lwd = 2)
    text(left + 0.25 * x_per_inch, key_y, key_labels[i], adj = c(0, 0.5), cex = 0.8)
    left <- left + key_widths[i] * x_per_inch
  }
  text(center, usr[4] + 0.3 * y_per_inch,
       "Sequence percent identity after T# (left to right): out-in1, out-in2, in1-in2; rounded to nearest integer (%)",
       cex = 0.8)
}

prepare_selected_tree_annotations <- function(tree, leaves, selected) {
  labeled <- tree
  labeled$tip.label <- accession_tip_labels(tree, leaves)
  roles <- trio_tip_roles(selected)
  tip_idx <- match(leaves$tip[match(roles$tip, leaves$tip)], tree$tip.label)
  if (anyNA(tip_idx)) {
    stop("Selected trio tip(s) missing from the tree or leaves: ",
         paste(roles$tip[is.na(tip_idx)], collapse = ", "), call. = FALSE)
  }

  annotated <- labeled
  for (i in seq_len(nrow(roles))) {
    idx <- tip_idx[i]
    annotated$tip.label[idx] <- paste0(annotated$tip.label[idx], "_trio",
                                      roles$trio[i], "_", roles$role[i])
  }

  ntip <- ape::Ntip(labeled)
  n_trio <- nrow(selected)
  tip_text <- paste0("T", roles$trio)
  tip_text[roles$role == "in1"] <- paste0(tip_text[roles$role == "in1"], ".1")
  tip_text[roles$role == "in2"] <- paste0(tip_text[roles$role == "in2"], ".2")
  out_members <- roles$role == "out"
  out_trios <- roles$trio[out_members]
  # Identities are nonnegative percentages; round .5 upward rather than to even.
  tip_text[out_members] <- sprintf(
    "T%d %.0f%% %.0f%% %.0f%%", out_trios,
    floor(as.numeric(selected$idt_12[out_trios]) + 0.5),
    floor(as.numeric(selected$idt_13[out_trios]) + 0.5),
    floor(as.numeric(selected$idt_23[out_trios]) + 0.5)
  )
  list(labeled = labeled, annotated = annotated, roles = roles, tip_idx = tip_idx,
       ntip = ntip, n_trio = n_trio, tip_text = tip_text)
}

write_selected_trio_tree <- function(tree, leaves, selected, out_dir, execution_info = NULL) {
  annotations <- prepare_selected_tree_annotations(tree, leaves, selected)
  labeled <- annotations$labeled; roles <- annotations$roles
  tip_idx <- annotations$tip_idx; tip_text <- annotations$tip_text
  ntip <- annotations$ntip; n_trio <- annotations$n_trio
  files <- stage1_tree_paths(out_dir)
  ape::write.tree(annotations$annotated, file = files[["nwk"]])
  layout <- selected_tree_pdf_layout(labeled, tip_idx, tip_text, n_trio, execution_info)
  open_tree_pdf(files[["pdf"]], layout$size)
  on.exit(grDevices::dev.off(), add = TRUE)
  margins <- par("mai")
  margins[3] <- margins[3] + layout$header_extra
  par(mai = margins)
  colors <- grDevices::hcl.colors(n_trio, "Dark 3")
  tip_colors <- ifelse(seq_len(ntip) %in% tip_idx, "black", "grey55")

  # Reserve physical space for each full annotation before the accession label.
  # Convert inches to tree units so short branches do not crowd the annotations.
  tree_depth <- if (is.null(tree$edge.length) || anyNA(tree$edge.length)) {
    max(ape::node.depth(tree)) - 1
  } else {
    max(ape::node.depth.edgelength(tree))
  }
  tree_depth <- max(tree_depth, .Machine$double.eps)
  annotation_width <- layout$annotation_width
  label_width <- layout$label_width
  units_per_inch <- tree_depth / max(par("pin")[[1]] - annotation_width - label_width, 1)
  label_offset <- annotation_width * units_per_inch
  ape::plot.phylo(
    labeled, cex = 0.7, underscore = TRUE, tip.color = tip_colors,
    label.offset = label_offset,
    x.lim = c(0, tree_depth + label_offset + label_width * units_per_inch),
    # Keep a gap below the legend, including for trees with many selected trios.
    y.lim = c(1, ntip + n_trio + 6), xaxs = "i"
  )

  # Repeated tips use successive positions, preserving the selected-trio order.
  tip_offsets <- numeric(ntip)
  for (k in seq_len(n_trio)) {
    members <- which(roles$trio == k)
    for (i in members) {
      idx <- tip_idx[i]
      offset <- (0.1 + tip_offsets[idx]) * units_per_inch
      ape::tiplabels(tip = idx, pch = 8, col = colors[k], cex = 1.6, lwd = 2,
                     offset = offset)
      if (roles$role[i] == "out") {
        ape::tiplabels(tip = idx, pch = 1, col = colors[k], cex = 1.6, lwd = 2,
                       offset = offset)
      }
      ape::tiplabels(text = tip_text[i], tip = idx, frame = "none",
                     col = colors[k], cex = 0.6, adj = c(0, 0.5),
                     offset = offset + 0.1 * units_per_inch)
      tip_offsets[idx] <- tip_offsets[idx] + layout$marker_steps[i]
    }
  }
  if (n_trio) {
    legend("topleft", legend = sprintf("T%d: %s (out) / %s + %s", seq_len(n_trio),
                                      selected$out_species, selected$in1_species,
                                      selected$in2_species),
           col = colors, pch = 8, bty = "n", cex = 0.7)
  }
  draw_selected_tree_header(n_trio, layout)
  invisible(files)
}
