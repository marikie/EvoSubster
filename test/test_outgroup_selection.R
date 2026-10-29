#!/usr/bin/env Rscript

# Tests for the outgroup-selection search in src/select/trio_selection.R:
#   - is_two_vs_one() genus pre-filter
#   - candidate_outgroups() ordering (nearest -> farthest)
#   - select_trios() early-stop: skip a too-close outgroup, pick the farther one,
#     stop at the first pass, and emit exactly one trio per ingroup pair.
#
# last-train is mocked (deterministic identities), so this runs with no network,
# no download, and no LAST -- only the real trio_filter.py judges each candidate.
#
# Run from the repo root:  Rscript test/test_outgroup_selection.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "src", "select", "trio_selection.R"))
# SCRIPT_DIR was resolved to test/ at source time; point TRIO_FILTER at the real script.
TRIO_FILTER <<- normalizePath(file.path(root, "src", "select", "trio_filter.py"))

suppressPackageStartupMessages(library(ape))

fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) {
    cat("ok:", msg, "\n")
  } else {
    cat("FAIL:", msg, "\n")
    fail <<- fail + 1L
  }
}

# --- is_two_vs_one ----------------------------------------------------------
check(is_two_vs_one("Homo", "Pan", "Pan") == FALSE,
      "two_vs_one: outgroup differs, ingroups share genus -> FALSE (pattern1)")
check(is_two_vs_one("Pan", "Pan", "Homo") == TRUE,
      "two_vs_one: outgroup congeneric with exactly one ingroup -> TRUE")
check(is_two_vs_one("Pan", "Pan", "Pan") == FALSE,
      "two_vs_one: all same genus -> FALSE (pattern3)")
check(is_two_vs_one("Gorilla", "Pan", "Homo") == FALSE,
      "two_vs_one: all different genera -> FALSE (pattern4)")

# --- candidate_outgroups ordering -------------------------------------------
tree <- read.tree(text = "(((A,B),C),D);")
leaves <- data.frame(
  tip = c("A", "B", "C", "D"),
  accession = c("a", "b", "c", "d"),
  species = c("Genus_a", "Genus_b", "Genus_c", "Genus_d"),
  genus = c("Genus", "Genus", "Genus", "Genus"),
  short_name = c("A", "B", "C", "D"),
  contig_n50 = c(1, 1, 1, 1),
  rel_contig_n50 = c(1, 1, 1, 1),
  stringsAsFactors = FALSE
)
check(identical(candidate_outgroups(tree, "A", "B", leaves), c("C", "D")),
      "candidate order for (A,B) is C then D (nearest -> farthest)")
check(length(candidate_outgroups(tree, "A", "D", leaves)) == 0,
      "pair whose MRCA is the root has no outgroup")

# --- ingroup_pairs: closest-first disjoint matching -------------------------
pair_keys <- function(pairs) sort(vapply(pairs, function(p) {
  s <- sort(p); paste(s[1], s[2], sep = "|")
}, character(1)))

# (((A,B),C),(D,E)): the two cherries (A,B) and (D,E) are the closest couples;
# both get matched, C is left unpaired (its only partners A/B are already used).
tree2 <- read.tree(text = "(((A,B),C),(D,E));")
leaves2 <- data.frame(
  tip = c("A", "B", "C", "D", "E"),
  accession = c("a", "b", "c", "d", "e"),
  species = c("Genus_a", "Genus_b", "Genus_c", "Genus_d", "Genus_e"),
  genus = c("Genus", "Genus", "Genus", "Genus", "Genus"),
  short_name = c("A", "B", "C", "D", "E"),
  contig_n50 = c(1, 1, 1, 1, 1),
  rel_contig_n50 = c(1, 1, 1, 1, 1),
  stringsAsFactors = FALSE
)
mpairs <- ingroup_pairs(tree2, leaves2, "matching")
check(identical(pair_keys(mpairs), c("A|B", "D|E")),
      "matching on (((A,B),C),(D,E)) yields the two cherries {A,B},{D,E}")
check({ tips_used <- unlist(mpairs); length(tips_used) == length(unique(tips_used)) },
      "matching uses each species at most once (disjoint)")
check(all(vapply(mpairs, function(p)
        length(candidate_outgroups(tree2, p[1], p[2], leaves2)) > 0, logical(1))),
      "every matched couple has an available outgroup (MRCA != root)")
check(length(setdiff(tree2$tip.label, unlist(mpairs))) == 1,
      "matching leaves exactly the stranded species (C) unpaired")

# "all" mode is the exhaustive C(n,2) enumeration (opt-out path).
check(length(ingroup_pairs(tree2, leaves2, "all")) == choose(5, 2),
      "all-pairing enumerates every tip pair")

# Regression: closeness must be topological, not branch-length.  Here the true cherry (A,B)
# has long terminal branches (patristic distance 6) while the NON-sister pair (A,C) has a
# SMALLER patristic distance (3.6); ranking by branch length would wrongly pair (A,C) and
# strand B.  Topological ranking (a cherry = 2 edges) must still pick the sister pair (A,B).
tree3 <- read.tree(text = "(((A:3,B:3):0.1,C:0.5):1,D:5);")
leaves3 <- data.frame(
  tip = c("A", "B", "C", "D"),
  accession = c("a", "b", "c", "d"),
  species = c("Genus_a", "Genus_b", "Genus_c", "Genus_d"),
  genus = c("Genus", "Genus", "Genus", "Genus"),
  short_name = c("A", "B", "C", "D"),
  contig_n50 = c(1, 1, 1, 1),
  rel_contig_n50 = c(1, 1, 1, 1),
  stringsAsFactors = FALSE
)
check(identical(pair_keys(ingroup_pairs(tree3, leaves3, "matching")), "A|B"),
      "topological matching picks the true cherry (A,B), not the shorter-patristic (A,C)")

# --- select_trios early-stop ------------------------------------------------
# For (A,B): nearest outgroup C is TOO CLOSE (idt_1x > idt_23, ordering fails);
# the farther outgroup D passes. Expect D chosen after 2 candidates tried.
idt_tbl <- c(
  "a|b" = 95,                 # ingroup pair
  "a|c" = 96, "b|c" = 96,     # C too close
  "a|d" = 90, "b|d" = 90      # D external -> passes the thesis rule
)
key <- function(x, y) if (x <= y) paste(x, y, sep = "|") else paste(y, x, sep = "|")
mock_get_identity <- function(acc_a, acc_b) unname(idt_tbl[key(acc_a, acc_b)])
fetch <- list(
  get_identity = mock_get_identity,
  counters = local({ e <- new.env(); e$trains <- 0L; e$downloads <- 0L; e })
)

opts <- modifyList(defaults, list(idt_threshold = 80, max_outgroup_tries = 5,
                                  out_dir = tempfile("ts_test_")))
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

sel <- select_trios(tree, leaves, opts, fetch = fetch)
ab <- sel[sel$in1_tip == "A" & sel$in2_tip == "B", , drop = FALSE]
check(nrow(ab) == 1, "(A,B) yields exactly one trio")
check(nrow(ab) == 1 && ab$out_tip == "D",
      "(A,B) outgroup is D (nearest C skipped as too-close)")
check(nrow(ab) == 1 && ab$outgroup_tries == "2",
      "(A,B) tried 2 candidates before the first pass")

# --- threshold options and their boundary behavior -------------------------
parsed <- parse_args(c("--tree", "example.nwk"))
check(parsed$ingroup_idt_threshold == 80 && parsed$outgroup_idt_threshold == 80 &&
        parsed$min_idt_gap_between_out_in == 0, "identity defaults are 80/80/0")
for (args in list(
  c("--idt-threshold", "99", "--ingroup-idt-threshold", "94", "--outgroup-idt-threshold", "89"),
  c("--ingroup-idt-threshold", "94", "--outgroup-idt-threshold", "89", "--idt-threshold", "99")
)) {
  parsed <- parse_args(c("--tree", "example.nwk", args))
  check(parsed$ingroup_idt_threshold == 94 && parsed$outgroup_idt_threshold == 89,
        "individual thresholds override common threshold regardless of argument order")
}
parsed <- parse_args(c("--tree", "example.nwk", "--idt-threshold", "85",
                      "--ingroup-idt-threshold", "94"))
check(parsed$ingroup_idt_threshold == 94 && parsed$outgroup_idt_threshold == 85,
      "common threshold remains the fallback for an unspecified individual threshold")
for (option in c("--idt-threshold", "--ingroup-idt-threshold", "--outgroup-idt-threshold",
                 "--min-idt-gap-between-out-in")) {
  for (value in c("-1", "100.1", "NaN", "Inf", "junk")) {
    error <- tryCatch(suppressWarnings(parse_args(c("--tree", "example.nwk", option, value))),
                      error = identity)
    check(inherits(error, "error") && grepl(option, conditionMessage(error), fixed = TRUE),
          paste("reject invalid identity option", option, value))
  }
  error <- tryCatch(parse_args(c("--tree", "example.nwk", option)), error = identity)
  check(inherits(error, "error"), paste("reject missing value for", option))
}

# C now has correct ordering, but one of its gaps is too small. D must be tried.
idt_tbl[c("a|c", "b|c")] <- c(92, 94)
gap_opts <- modifyList(opts, list(min_idt_gap_between_out_in = 2))
sel <- select_trios(tree, leaves, gap_opts, fetch = fetch)
check(nrow(sel) == 1 && sel$out_tip == "D" && sel$outgroup_tries == "2",
      "a candidate with only one sufficient gap is skipped in favor of the next")
check(nrow(sel) == 1 && as.numeric(sel$idt_gap_12) == 5 && as.numeric(sel$idt_gap_13) == 5 &&
        as.numeric(sel$min_idt_gap_between_out_in_value) == 2,
      "selected trio records both gaps and the applied gap threshold")
sel <- select_trios(tree, leaves, modifyList(gap_opts, list(min_idt_gap_between_out_in = 5)), fetch)
check(nrow(sel) == 0, "gap equal to configured threshold fails")

# An absolute outgroup failure closes the search window before D is trained.
idt_tbl[c("a|c", "b|c")] <- c(89, 89)
seen <- character()
counting_fetch <- fetch
counting_fetch$get_identity <- function(a, b) {
  seen <<- c(seen, key(a, b))
  mock_get_identity(a, b)
}
sel <- select_trios(tree, leaves, modifyList(opts, list(outgroup_idt_threshold = 89)), counting_fetch)
check(nrow(sel) == 0 && !any(c("a|d", "b|d") %in% seen),
      "outgroup threshold equality stops the search before farther candidates")

idt_tbl[c("a|b", "a|c", "b|c")] <- c(79, 75, 75)
sel <- select_trios(tree, leaves,
                    modifyList(opts, list(ingroup_idt_threshold = 78, outgroup_idt_threshold = 74)), fetch)
check(nrow(sel) == 1 && sel$out_tip == "C", "independent low thresholds reach the Python verdict")
sel <- select_trios(tree, leaves,
                    modifyList(opts, list(ingroup_idt_threshold = 79, outgroup_idt_threshold = 74)), fetch)
check(nrow(sel) == 0, "ingroup threshold equality fails the prefilter")

exempt_leaves <- leaves
exempt_leaves$genus[exempt_leaves$tip == "C"] <- "Other"
sel <- select_trios(tree, exempt_leaves,
                    modifyList(opts, list(ingroup_idt_threshold = 78, outgroup_idt_threshold = 74,
                                          min_idt_gap_between_out_in = 20)), fetch)
check(nrow(sel) == 1 && sel$out_tip == "C" && sel$idt_gap_exempt == "TRUE" &&
        sel$idt_gap_condition == "FALSE", "genus exemption survives nonzero gap selection")

if (fail > 0L) {
  cat("\n", fail, "test(s) FAILED\n", sep = "")
  quit(status = 1)
}
cat("\nAll outgroup-selection tests passed.\n")
