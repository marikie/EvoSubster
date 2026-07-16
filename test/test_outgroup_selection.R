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
  stringsAsFactors = FALSE
)
check(identical(pair_keys(ingroup_pairs(tree3, leaves3, "matching")), "A|B"),
      "topological matching picks the true cherry (A,B), not the shorter-patristic (A,C)")

# --- select_trios early-stop ------------------------------------------------
# For (A,B): nearest outgroup C is TOO CLOSE (idt_1x > idt_23, ordering fails);
# the farther outgroup D passes. Expect D chosen after 2 candidates tried.
idt_tbl <- c(
  "A|B" = 95,                 # ingroup pair
  "A|C" = 96, "B|C" = 96,     # C too close
  "A|D" = 90, "B|D" = 90      # D external -> passes the thesis rule
)
key <- function(x, y) if (x <= y) paste(x, y, sep = "|") else paste(y, x, sep = "|")
mock_get_identity <- function(sa, aa, sb, ab) unname(idt_tbl[key(sa, sb)])
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

if (fail > 0L) {
  cat("\n", fail, "test(s) FAILED\n", sep = "")
  quit(status = 1)
}
cat("\nAll outgroup-selection tests passed.\n")
