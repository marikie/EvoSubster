#!/usr/bin/env Rscript

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "src", "select", "trio_selection.R"))

entry <- function(info, key) info$settings[info$settings$key == key, , drop = FALSE]
args <- c("--tree", "input.nwk")
info <- record_selection_execution(parse_args(args), args, "src/select/trio_selection.R")
stopifnot(
  identical(info$command, "Rscript src/select/trio_selection.R --tree input.nwk"),
  setequal(info$settings$key, names(defaults)),
  entry(info, "ingroup_idt_threshold")$value == "80",
  entry(info, "ingroup_idt_threshold")$source == "default",
  entry(info, "outgroup_idt_threshold")$value == "80",
  entry(info, "min_idt_gap_between_out_in")$value == "0",
  entry(info, "max_outgroup_tries")$value == "5",
  entry(info, "ingroup_pairing")$value == "matching",
  entry(info, "min_rel_contig_n50")$value == "0 (off)",
  entry(info, "stage0_top_k")$value == "3",
  entry(info, "assembly_qc")$value == "none",
  entry(info, "allow_qc_override")$value == "off",
  entry(info, "dry_run")$value == "off",
  entry(info, "keep_unused_species_data")$value == "off",
  entry(info, "date")$value == format(Sys.Date(), "%Y%m%d"),
  entry(info, "train_cache_dir")$value == "./results/trio_selection/train_cache",
  entry(info, "train_cache_dir")$source == "derived from --out-dir"
)

args <- c("--tree", "tree with 'quote' $(printf BAD).nwk",
          "--ingroup-idt-threshold", "80", "--idt-threshold", "85.5",
          "--min-idt-gap-between-out-in", "1", "--threads", "8",
          "--date", "20260102", "--allow-qc-override", "--keep-unused-species-data",
          "--out-dir", "output with spaces")
info <- record_selection_execution(parse_args(args), args, "path with spaces/trio_selection.R")
stopifnot(
  entry(info, "ingroup_idt_threshold")$value == "80",
  entry(info, "ingroup_idt_threshold")$source == "specified",
  entry(info, "outgroup_idt_threshold")$value == "85.5",
  entry(info, "outgroup_idt_threshold")$source == "from --idt-threshold",
  entry(info, "idt_threshold")$source == "specified",
  entry(info, "threads")$source == "specified",
  entry(info, "allow_qc_override")$value == "on",
  entry(info, "keep_unused_species_data")$value == "on",
  entry(info, "date")$value == "20260102",
  entry(info, "train_cache_dir")$value == "output with spaces/train_cache"
)
roundtrip <- system2("sh", c("-c", shQuote(paste("printf '%s\\n'", info$command))), stdout = TRUE)
stopifnot(identical(roundtrip, c("Rscript", "path with spaces/trio_selection.R", args)))
for (options in list(c("--idt-threshold", "85", "--outgroup-idt-threshold", "80"),
                     c("--outgroup-idt-threshold", "80", "--idt-threshold", "85"))) {
  argv <- c("--tree", "input.nwk", options, "--train-cache-dir", "cache-dir")
  result <- record_selection_execution(parse_args(argv), argv, "trio_selection.R")
  stopifnot(entry(result, "outgroup_idt_threshold")$source == "specified",
            entry(result, "ingroup_idt_threshold")$source == "from --idt-threshold",
            entry(result, "train_cache_dir")$source == "specified",
            entry(result, "train_cache_dir")$value == "cache-dir")
}
stopifnot(identical(execution_info_lines(NULL),
                    c("Command: not recorded", "Effective settings: not recorded")))
cat("ok: all effective values, option provenance and shell quoting\n")

if (!requireNamespace("ape", quietly = TRUE)) stop("PDF tests require ape")
if (!nzchar(Sys.which("pdftotext"))) stop("PDF tests require pdftotext")
work <- tempfile("selection-execution-")
dir.create(work)
tree <- ape::read.tree(text = "((A:1,B:1):1,C:2);")
leaves <- data.frame(tip = c("A", "B", "C"), accession = c("a", "b", "c"))
selected <- data.frame(out_tip = "C", in1_tip = "A", in2_tip = "B",
                       out_species = "Other c", in1_species = "Genus a", in2_species = "Genus b",
                       idt_12 = 97.1778, idt_13 = 97.5568, idt_23 = 97.6227)

pdf_text <- function(path) paste(system2("pdftotext", c(shQuote(path), "-"), stdout = TRUE), collapse = "\n")
files <- write_selected_trio_tree(tree, leaves, selected, work, execution_info = info)
contents <- pdf_text(files[["pdf"]])
for (phrase in c("Command:", "Rscript", "--min-idt-gap-between-out-in", "Effective selection settings:",
                 "Ingroup identity > 80% [specified]", "Outgroup-ingroup identity > 85.5%",
                 "Identity gap > 1 percentage point", "[from --idt-threshold]",
                 "Gap exemption:", "Runtime:", "Paths:", "Threads: 8 [specified]", "T1 97% 98% 98%")) {
  stopifnot(grepl(phrase, contents, fixed = TRUE))
}
stopifnot(!grepl("idt-threshold:", contents, fixed = TRUE))
short_pdf <- file.path(work, "short.pdf")
file.copy(files[["pdf"]], short_pdf)

# A single long path must wrap without truncation, not shrink the tree area.
long_path <- paste0("/", paste(rep("long_path_segment", 70), collapse = "/"), "/input.nwk")
argv <- c("--tree", long_path, "--ingroup-idt-threshold", "90")
long_info <- record_selection_execution(parse_args(argv), argv, "src/select/trio_selection.R")
files <- write_selected_trio_tree(tree, leaves, selected, work, execution_info = long_info)
contents <- pdf_text(files[["pdf"]])
stopifnot(grepl(long_path, gsub("[[:space:]]", "", contents), fixed = TRUE))
long_pdf <- file.path(work, "long.pdf")
file.copy(files[["pdf"]], long_pdf)

pdf_geometry <- function(path) {
  xml <- system2("pdftotext", c("-bbox", shQuote(path), "-"), stdout = TRUE)
  page <- regmatches(xml, regexec('<page width="([0-9.]+)" height="([0-9.]+)">', xml))
  page <- Filter(function(x) length(x) == 3, page)
  stopifnot(length(page) == 1)
  words <- regmatches(xml, regexec(
    '<word xMin="([0-9.]+)" yMin="([0-9.]+)" xMax="([0-9.]+)" yMax="([0-9.]+)">(.*)</word>', xml))
  words <- Filter(function(x) length(x) == 6, words)
  bounds <- t(vapply(words, function(x) as.numeric(x[2:5]), numeric(4)))
  stopifnot(all(bounds[, 1:2] >= 0), all(bounds[, 3] <= as.numeric(page[[1]][2])),
            all(bounds[, 4] <= as.numeric(page[[1]][3])))
  list(height = as.numeric(page[[1]][3]), bounds = bounds,
       text = vapply(words, function(x) x[6], character(1)))
}
short_geometry <- pdf_geometry(short_pdf)
long_geometry <- pdf_geometry(long_pdf)
stopifnot(long_geometry$height > short_geometry$height)
for (geometry in list(short_geometry, long_geometry)) {
  y <- function(word, edge) geometry$bounds[geometry$text == word, edge]
  stopifnot(min(y("Stage", 4)) < min(y("Command:", 2)),
            max(y("Assembly", 4)) < min(y("Outgroup", 2)),
            max(y("Sequence", 4)) < min(y("A_a", 2)))
}
files <- write_selected_trio_tree(tree, leaves, selected, work)
stopifnot(grepl("Command: not recorded", pdf_text(files[["pdf"]]), fixed = TRUE))

# Check that main passes the actual invocation and parsed options to the PDF.
local({
  original_prune <- prune_to_best_assembly
  original_fetch <- make_fetchers
  original_select <- select_trios
  on.exit({
    prune_to_best_assembly <<- original_prune
    make_fetchers <<- original_fetch
    select_trios <<- original_select
  })
  input <- file.path(work, "main input.nwk")
  ape::write.tree(tree, input)
  main_leaves <- leaves
  main_leaves$species <- c("Genus a", "Genus b", "Other c")
  main_leaves$rel_contig_n50 <- 0
  prune_to_best_assembly <<- function(...) list(tree = tree, leaves = main_leaves)
  make_fetchers <<- function(...) list()
  select_trios <<- function(...) selected
  main(c("--tree", input, "--out-dir", work, "--ingroup-idt-threshold", "90",
         "--keep-unused-species-data"), script_path = "src/select/trio_selection.R")
  contents <- pdf_text(file.path(work, "selected_trios.pdf"))
  stopifnot(grepl("Ingroup identity > 90% [specified]", contents, fixed = TRUE),
            grepl("Outgroup-ingroup identity > 80% [default]", contents, fixed = TRUE),
            grepl("Rscript src/select/trio_selection.R", contents, fixed = TRUE))
})

if (nzchar(Sys.getenv("SELECTION_PREVIEW_DIR"))) {
  dir.create(Sys.getenv("SELECTION_PREVIEW_DIR"), recursive = TRUE, showWarnings = FALSE)
  file.copy(c(short_pdf, long_pdf), Sys.getenv("SELECTION_PREVIEW_DIR"), overwrite = TRUE)
}
unlink(work, recursive = TRUE)
cat("ok: command and settings reach the final PDF, including long paths and legacy callers\n")
