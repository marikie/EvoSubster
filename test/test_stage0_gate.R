#!/usr/bin/env Rscript

# Tests for Stage 0 assembly selection in src/select/trio_selection.R:
#   select_best_assemblies() -- current-only filter, relative contig-N50 gate, and
#   reference-first best-per-species pick.  Runs with no network (synthetic metadata).
#
# Run from the repo root:  Rscript test/test_stage0_gate.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "src", "select", "assembly_selection.R"))

fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) cat("ok:", msg, "\n") else { cat("FAIL:", msg, "\n"); fail <<- fail + 1L }
}

# One row per assembly.  total_ungapped_length = 1e8 throughout, so
# rel_contig_n50 = contig_n50 / 1e8.
meta <- data.frame(
  accession = c("A_ref", "A_big", "B_big", "B_small", "G_ref", "D_ref_supp", "D_alt",
                "E_annot", "E_big"),
  organism_name = c("Genus alpha", "Genus alpha", "Genus beta", "Genus beta",
                    "Genus gamma", "Genus delta", "Genus delta",
                    "Genus epsilon", "Genus epsilon"),
  species = c("Genus alpha", "Genus alpha", "Genus beta", "Genus beta",
              "Genus gamma", "Genus delta", "Genus delta",
              "Genus epsilon", "Genus epsilon"),
  genus = "Genus",
  refseq_category = c("reference genome", "", "", "", "reference genome", "reference genome", "",
                      "", ""),
  assembly_level = c("Chromosome", "Contig", "Contig", "Contig", "Chromosome", "Chromosome", "Contig",
                     "Chromosome", "Contig"),
  contig_n50 = c(1e6, 5e6, 2e6, 1e5, 1e4, 1e6, 2e6,
                 1e6, 4e6),
  total_ungapped_length = 1e8,
  total_number_of_chromosomes = 10,
  has_annotation = c("true", "false", "false", "false", "true", "true", "false",
                     "true", "false"),
  assembly_status = c("current", "current", "current", "current", "current", "suppressed", "current",
                      "current", "current"),
  stringsAsFactors = FALSE
)

pick <- function(best, sp) best$accession[best$species == sp]

# Gate off: every species with >=1 current assembly is kept.
b0 <- select_best_assemblies(meta, 0)
check(setequal(b0$species, c("Genus alpha", "Genus beta", "Genus gamma", "Genus delta",
                             "Genus epsilon")),
      "gate off: all five species kept")
check(identical(pick(b0, "Genus alpha"), "A_ref"),
      "reference-first: alpha picks the reference (A_ref) over the more contiguous A_big")
check(identical(pick(b0, "Genus beta"), "B_big"),
      "no reference: beta falls back to the most contiguous (B_big)")
check(identical(pick(b0, "Genus delta"), "D_alt"),
      "suppressed excluded: delta picks the current D_alt, not the suppressed reference")
check(identical(pick(b0, "Genus epsilon"), "E_annot"),
      "general fallback: chromosome-level E_annot beats the more fragmented contig-level E_big")

# Gate on at rel N50 >= 0.005 (= 5e5 bp at 1e8): gamma's only assembly (1e4 -> 1e-4) fails.
b1 <- select_best_assemblies(meta, 0.005)
check(!("Genus gamma" %in% b1$species),
      "gate on: gamma dropped (its only assembly is below the relative floor)")
check(setequal(b1$species, c("Genus alpha", "Genus beta", "Genus delta", "Genus epsilon")),
      "gate on: the other four species survive")
check(identical(pick(b1, "Genus beta"), "B_big"),
      "gate on: B_small (rel 0.001) filtered, beta still picks B_big")

if (fail > 0L) { cat("\n", fail, " test(s) FAILED\n", sep = ""); quit(status = 1) }
cat("\nAll Stage 0 gate tests passed.\n")
