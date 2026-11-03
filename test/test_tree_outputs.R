#!/usr/bin/env Rscript

# Unit tests for trio-selection Newick and PDF outputs using ape.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))

if (!requireNamespace("ape", quietly = TRUE)) {
  cat("SKIP: tree-output tests require ape, which is not installed.\n")
  quit(status = 0)
}

source(file.path(root, "src", "select", "trio_selection.R"))

fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) {
    cat("ok:", msg, "\n")
  } else {
    cat("FAIL:", msg, "\n")
    fail <<- fail + 1L
  }
}

tree <- ape::read.tree(
  text = "((A_sp:1,B_sp:1):1,((C_sp:1,D_sp:1):1,(E_sp:1,F_sp:1):1):1);"
)
leaves <- data.frame(
  tip = c("A_sp", "B_sp", "C_sp", "D_sp", "E_sp", "F_sp"),
  species = c("A sp", "B sp", "C sp", "D sp", "E sp", "F sp"),
  accession = sprintf("GCA_%09d.1", 1:6),
  short_name = c("Asp1", "Bsp2", "Csp3", "Dsp4", "Esp5", "Fsp6"),
  stringsAsFactors = FALSE
)
# A different metadata order catches accidental positional matching.
leaves <- leaves[c(4, 2, 6, 1, 5, 3), , drop = FALSE]
base_labels <- paste0(tree$tip.label, "_", sprintf("GCA_%09d.1", 1:6))

check(
  identical(accession_tip_labels(tree, leaves), base_labels),
  "accession labels follow tree tip order when metadata rows are reordered"
)
missing_tip_error <- try(
  accession_tip_labels(tree, leaves[leaves$tip != "C_sp", , drop = FALSE]),
  silent = TRUE
)
check(
  inherits(missing_tip_error, "try-error") &&
    grepl("C_sp", as.character(missing_tip_error), fixed = TRUE),
  "accession labels report the tree tip missing from leaves"
)
check(
  identical(tree_pdf_size(6L), c(width = 11, height = 7)) &&
    identical(tree_pdf_size(100L), c(width = 11, height = 20)),
  "PDF dimensions retain the minimum height and grow with tree size"
)

stage0_dir <- tempfile("trio-tree-stage0-", tmpdir = tempdir())
dir.create(stage0_dir)
device_before <- grDevices::dev.cur()
stage0_result <- withVisible(write_stage0_tree(tree, leaves, stage0_dir))
stage0_files <- stage0_result$value
check(
  !stage0_result$visible &&
    identical(stage0_files, c(
      nwk = file.path(stage0_dir, "pruned_tree.nwk"),
      pdf = file.path(stage0_dir, "pruned_tree.pdf")
    )),
  "Stage 0 returns the named Newick and PDF paths invisibly"
)
stage0_tree <- ape::read.tree(stage0_files[["nwk"]])
check(
  ape::Ntip(stage0_tree) == ape::Ntip(tree) &&
    identical(stage0_tree$tip.label, accession_tip_labels(tree, leaves)),
  "Stage 0 Newick preserves every tip with its accession label"
)
check(
  file.exists(stage0_files[["pdf"]]) &&
    file.info(stage0_files[["pdf"]])$size > 0,
  "Stage 0 PDF exists and is nonempty"
)
check(
  identical(readChar(stage0_files[["pdf"]], nchars = 4L, useBytes = TRUE), "%PDF"),
  "Stage 0 PDF starts with the PDF signature"
)
check(
  identical(grDevices::dev.cur(), device_before),
  "Stage 0 closes its PDF device"
)

selected <- data.frame(
  out_tip = c("A_sp", "D_sp"),
  in1_tip = c("B_sp", "E_sp"),
  in2_tip = c("C_sp", "F_sp"),
  out_acc = sprintf("GCA_%09d.1", c(1, 4)),
  in1_acc = sprintf("GCA_%09d.1", c(2, 5)),
  in2_acc = sprintf("GCA_%09d.1", c(3, 6)),
  out_species = c("A sp", "D sp"),
  in1_species = c("B sp", "E sp"),
  in2_species = c("C sp", "F sp"),
  out_short = c("Asp1", "Dsp4"),
  in1_short = c("Bsp2", "Esp5"),
  in2_short = c("Csp3", "Fsp6"),
  genus_1 = c("A", "D"),
  genus_2 = c("B", "E"),
  genus_3 = c("C", "F"),
  idt_12 = c(85, 86),
  idt_13 = c(84, 85),
  idt_23 = c(95, 96),
  outgroup_tries = c(1L, 2L),
  stringsAsFactors = FALSE
)
roles <- trio_tip_roles(selected)
check(
  is.data.frame(roles) && identical(names(roles), c("tip", "trio", "role")) &&
    nrow(roles) == 6L && is.integer(roles$trio) &&
    identical(sort(unique(roles$trio)), 1:2),
  "two selected rows expand into six role rows with integer trio numbers"
)
for (k in seq_len(nrow(selected))) {
  trio_roles <- roles[roles$trio == k, , drop = FALSE]
  check(
    nrow(trio_roles) == 3L &&
      setequal(trio_roles$role, c("out", "in1", "in2")) &&
      identical(
        trio_roles$tip[match(c("out", "in1", "in2"), trio_roles$role)],
        c(selected$out_tip[k], selected$in1_tip[k], selected$in2_tip[k])
      ),
    sprintf("trio %d roles match selected row %d", k, k)
  )
}
for (empty_selected in list(selected[FALSE, , drop = FALSE], data.frame())) {
  empty_roles <- trio_tip_roles(empty_selected)
  check(
    is.data.frame(empty_roles) && nrow(empty_roles) == 0L &&
      identical(names(empty_roles), c("tip", "trio", "role")) &&
      is.integer(empty_roles$trio),
    "empty selected tables produce the zero-row tip/trio/role schema"
  )
}

stage1_dir <- tempfile("trio-tree-stage1-", tmpdir = tempdir())
dir.create(stage1_dir)
stage1_result <- withVisible(write_selected_trio_tree(
  tree, leaves, selected[1, , drop = FALSE], stage1_dir
))
stage1_files <- stage1_result$value
check(
  !stage1_result$visible &&
    identical(stage1_files, c(
      nwk = file.path(stage1_dir, "selected_trios.nwk"),
      pdf = file.path(stage1_dir, "selected_trios.pdf")
    )),
  "Stage 1 returns the named Newick and PDF paths invisibly"
)
stage1_tree <- ape::read.tree(stage1_files[["nwk"]])
expected_labels <- base_labels
expected_labels[1:3] <- paste0(
  expected_labels[1:3], c("_trio1_out", "_trio1_in1", "_trio1_in2")
)
check(
  ape::Ntip(stage1_tree) == ape::Ntip(tree) &&
    identical(stage1_tree$tip.label, expected_labels),
  "Stage 1 Newick appends the correct trio and role suffix to each member"
)
check(
  identical(stage1_tree$tip.label[4:6], base_labels[4:6]) &&
    !any(grepl("_trio[0-9]+_(out|in1|in2)", stage1_tree$tip.label[4:6])),
  "tips outside every trio retain their base accession labels"
)
check(
  file.exists(stage1_files[["pdf"]]) &&
    file.info(stage1_files[["pdf"]])$size > 0,
  "Stage 1 PDF exists and is nonempty"
)
check(
  identical(grDevices::dev.cur(), device_before),
  "Stage 1 closes its PDF device"
)

if (nzchar(Sys.which("pdffonts"))) {
  pdf_paths <- c("Stage 0" = stage0_files[["pdf"]], "Stage 1" = stage1_files[["pdf"]])
  for (stage in names(pdf_paths)) {
    path <- pdf_paths[[stage]]
    rows <- system2("pdffonts", shQuote(path), stdout = TRUE)
    rows <- rows[-seq_len(2L)]
    font_flags <- regmatches(rows, regexec(
      "(yes|no)\\s+(yes|no)\\s+(yes|no)\\s+[0-9]+\\s+[0-9]+\\s*$",
      rows, perl = TRUE
    ))
    check(
      length(rows) >= 1L && all(vapply(font_flags, function(flags) {
        length(flags) == 4L && flags[[2]] == "yes"
      }, logical(1))),
      paste(stage, "PDF embeds every font")
    )
  }
} else {
  cat("skipping pdffonts embedding check\n")
}

stage1_stage0_files <- write_stage0_tree(tree, leaves, stage1_dir)
stage0_checksums <- tools::md5sum(stage1_stage0_files)
removed <- withVisible(remove_stage1_tree_outputs(stage1_dir))
check(
  !removed$visible && identical(removed$value, stage1_files),
  "Stage 1 cleanup returns the named Newick and PDF paths invisibly"
)
check(
  !any(file.exists(stage1_files)),
  "Stage 1 cleanup removes both previously written tree files"
)
check(
  all(file.exists(stage1_stage0_files)) &&
    identical(tools::md5sum(stage1_stage0_files), stage0_checksums),
  "Stage 1 cleanup leaves pruned_tree.nwk and pruned_tree.pdf unchanged"
)
missing_cleanup <- try(remove_stage1_tree_outputs(stage1_dir), silent = TRUE)
check(
  !inherits(missing_cleanup, "try-error") &&
    identical(missing_cleanup, stage1_files),
  "Stage 1 cleanup succeeds when both tree files are already absent"
)

shared <- selected
shared[2, c("out_tip", "out_acc", "out_species", "out_short", "genus_1")] <-
  shared[1, c("out_tip", "out_acc", "out_species", "out_short", "genus_1")]
shared_dir <- tempfile("trio-tree-shared-", tmpdir = tempdir())
dir.create(shared_dir)
shared_files <- write_selected_trio_tree(tree, leaves, shared, shared_dir)
shared_tree <- ape::read.tree(shared_files[["nwk"]])
shared_labels <- base_labels
shared_labels[c(1, 2, 3, 5, 6)] <- paste0(
  shared_labels[c(1, 2, 3, 5, 6)],
  c("_trio1_out_trio2_out", "_trio1_in1", "_trio1_in2", "_trio2_in1", "_trio2_in2")
)
check(
  identical(shared_tree$tip.label, shared_labels) &&
    endsWith(shared_tree$tip.label[1], "_trio1_out_trio2_out"),
  "a shared outgroup accumulates suffixes in selected trio order"
)
check(
  file.exists(shared_files[["pdf"]]) &&
    file.info(shared_files[["pdf"]])$size > 0 &&
    identical(grDevices::dev.cur(), device_before),
  "shared-outgroup PDF is written and its device is closed"
)
check(
  identical(tree$tip.label, c("A_sp", "B_sp", "C_sp", "D_sp", "E_sp", "F_sp")),
  "writing tree outputs leaves the original tree labels unchanged"
)

dry_dir <- file.path(tempdir(), "dryrun")
dir.create(dry_dir, recursive = TRUE)
tree_file <- file.path(dry_dir, "input.nwk")
ape::write.tree(tree, file = tree_file)
stale <- stage1_tree_paths(dry_dir)
for (p in stale) writeLines("stale", p)
local({
  saved_parse <- parse_args
  saved_prune <- prune_to_best_assembly
  on.exit({
    parse_args <<- saved_parse
    prune_to_best_assembly <<- saved_prune
  }, add = TRUE)
  dry_leaves <- leaves
  dry_leaves$rel_contig_n50 <- 0
  parse_args <<- function(argv) modifyList(defaults, list(
    tree = tree_file, out_dir = dry_dir, dry_run = TRUE
  ))
  prune_to_best_assembly <<- function(tree, ...) list(tree = tree, leaves = dry_leaves)
  main()
})
check(
  !any(file.exists(stale)),
  "--dry-run removes Stage 1 tree files left by a previous run"
)
dry_outputs <- file.path(dry_dir, c(
  "pruned_tree.nwk", "pruned_tree.pdf",
  "selected_assemblies.tsv", "candidate_outgroups.tsv"
))
check(
  all(file.exists(dry_outputs)) && all(file.info(dry_outputs)$size > 0),
  "--dry-run writes nonempty Stage 0 tree files and assembly and candidate TSVs"
)
check(
  grDevices::dev.cur() == 1,
  "--dry-run leaves no graphics device open"
)

if (fail > 0L) {
  cat("\n", fail, " test(s) FAILED\n", sep = "")
  quit(status = 1)
}

cat("\nAll tree-output tests passed.\n")
