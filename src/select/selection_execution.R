# Execution metadata and text layout for the final selection tree.

record_selection_execution <- function(opts, argv, script_path) {
  quote_arg <- function(value) {
    if (nzchar(value) && grepl("^[A-Za-z0-9_./:=,+%@-]+$", value)) value else shQuote(value, type = "sh")
  }
  command <- if (is.null(script_path) || !nzchar(script_path)) "not recorded" else {
    paste(vapply(c("Rscript", script_path, argv), quote_arg, character(1)), collapse = " ")
  }
  specified <- attr(opts, "specified_options")
  keys <- c("idt_threshold", "ingroup_idt_threshold", "outgroup_idt_threshold",
            "min_idt_gap_between_out_in", "ingroup_pairing", "max_outgroup_tries",
            "min_rel_contig_n50", "stage0_top_k", "allow_qc_override",
            "threads", "date", "dry_run", "keep_unused_species_data",
            "tree", "genome_dir", "out_dir", "train_cache_dir", "assembly_qc")
  labels <- c("Common identity threshold (fallback)", "Ingroup identity", "Outgroup-ingroup identity",
              "Identity gap", "Pairing", "Max outgroup tries", "Relative contig N50 minimum",
              "Stage 0 top-k", "QC override", "Threads", "Date", "Dry run",
              "Keep unused species data", "Tree", "Genome directory", "Output directory",
              "Train cache", "Assembly QC")
  sections <- c(rep("Effective selection settings", 9), rep("Runtime", 4), rep("Paths", 5))
  settings <- data.frame(key = keys, section = sections, label = labels,
                         value = "", source = "", stringsAsFactors = FALSE)
  for (i in seq_along(keys)) {
    key <- keys[i]
    value <- opts[[key]]
    origin <- if (is.null(specified)) "not recorded" else if (key %in% specified) "specified" else "default"
    if (key %in% c("ingroup_idt_threshold", "outgroup_idt_threshold") &&
        !key %in% specified && "idt_threshold" %in% specified) origin <- "from --idt-threshold"
    if (key == "train_cache_dir" && (is.null(value) || !nzchar(value))) {
      value <- file.path(opts$out_dir, "train_cache")
      origin <- "derived from --out-dir"
    }
    settings$value[i] <- if (is.null(value)) {
      "none"
    } else if (is.logical(value)) {
      if (value) "on" else "off"
    } else if (is.numeric(value)) {
      format(value, trim = TRUE, scientific = FALSE, digits = 15)
    } else as.character(value)
    if (key == "min_rel_contig_n50" && identical(as.numeric(value), 0)) settings$value[i] <- "0 (off)"
    settings$source[i] <- origin
  }
  list(command = command, settings = settings)
}

execution_info_lines <- function(execution_info) {
  if (is.null(execution_info)) return(c("Command: not recorded", "Effective settings: not recorded"))
  lines <- c("Command:", paste0("  ", execution_info$command))
  for (section in unique(execution_info$settings$section)) {
    lines <- c(lines, "", paste0(section, ":"))
    rows <- execution_info$settings[execution_info$settings$section == section, , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      key <- rows$key[i]
      value <- rows$value[i]
      label <- rows$label[i]
      text <- if (key %in% c("ingroup_idt_threshold", "outgroup_idt_threshold")) {
        paste0(label, " > ", value, "%")
      } else if (key == "min_idt_gap_between_out_in") {
        unit <- if (value == "1") "percentage point" else "percentage points"
        paste0(label, " > ", value, " ", unit, " (both outgroup-ingroup comparisons)")
      } else if (key == "idt_threshold") {
        paste0(label, ": ", value, "%")
      } else paste0(label, ": ", value)
      lines <- c(lines, paste0("  ", text, " [", rows$source[i], "]"))
    }
    if (section == "Effective selection settings") {
      lines <- c(lines, "  Gap exemption: ingroups share a genus distinct from the outgroup")
    }
  }
  lines
}

wrap_execution_lines <- function(lines, width, cex = 0.72) {
  # Font measurements, including hard breaks for a single oversized path token.
  wrapped <- character()
  for (line in lines) {
    line <- gsub("\n", "\\n", line, fixed = TRUE)
    while (nzchar(line) && strwidth(line, units = "inches", cex = cex) > width) {
      lo <- 1L; hi <- nchar(line)
      while (lo < hi) {
        mid <- ceiling((lo + hi) / 2)
        if (strwidth(substr(line, 1, mid), units = "inches", cex = cex) <= width) lo <- mid else hi <- mid - 1L
      }
      spaces <- gregexpr(" ", substr(line, 1, lo), fixed = TRUE)[[1]]
      spaces <- spaces[spaces > 2]
      cut <- if (length(spaces)) tail(spaces, 1) else lo
      wrapped <- c(wrapped, trimws(substr(line, 1, cut), which = "right"))
      line <- trimws(substring(line, cut + 1), which = "left")
    }
    wrapped <- c(wrapped, line)
  }
  wrapped
}
