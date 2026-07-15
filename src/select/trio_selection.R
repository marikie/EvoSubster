#!/usr/bin/env Rscript

# Select substitution-analysis trios from a Newick tree.
#
# Original script written by Ryota Ishii; rewritten to apply the thesis
# (section 2.3) rule -- measured sequence percent identity plus genus names --
# instead of tree distances.
#
# The tree cannot decide a trio on its own.  On an ultrametric tree the outgroup
# sits at an identical distance from both ingroup species, so the thesis'
# identity ordering test is satisfied by *every* trio the tree can produce and
# filters nothing.  The asymmetry that test looks for only exists in real aligned
# sequence, so the tree is used to enumerate and cheaply prescreen candidates,
# and last-train supplies the identities that actually decide them.
#
# Stage 0  prune the tree to one assembly per species   [NCBI metadata, no download]
# Stage 1  enumerate trios and prescreen on the tree    [free]
# Stage 2  download the surviving accessions            [dwl_organism.sh]
# Stage 3  last-train once per unique species pair      [last_train.sh, shared cache]
# Stage 4  apply the thesis rule                        [trio_filter.py]
#
# Usage:
#   Rscript src/select/trio_selection.R --tree Primates_tree.nwk [options]
#
# Depends on: ape, dplyr (R); python3, LAST, curl, jq (external).

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
})

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
})
SRC_DIR <- dirname(SCRIPT_DIR)

FETCH_METADATA <- file.path(SCRIPT_DIR, "fetch_assembly_metadata.py")
TRIO_FILTER <- file.path(SCRIPT_DIR, "trio_filter.py")
DWL_ORGANISM <- file.path(SRC_DIR, "dwl_organism.sh")
LAST_TRAIN <- file.path(SRC_DIR, "align", "last_train.sh")

IDENTITY_PREFIX <- "# substitution percent identity:"

# Best assembly per species: prefer the designated reference, then the most
# contiguous assembly.  A strict reference-only filter would drop whole species
# from the tree (11 of 37 primates), which is worse than keeping a good
# non-reference assembly.
ASSEMBLY_LEVEL_RANK <- c(
  "Complete Genome" = 1,
  "Chromosome" = 2,
  "Scaffold" = 3,
  "Contig" = 4
)

# --- options ----------------------------------------------------------------

defaults <- list(
  tree = NULL,
  idt_threshold = 80,
  prescreen_margin = 10,
  top_k = 3,
  ingroup_min = 0,
  min_contig_n50 = 0,
  genome_dir = "./genomes",
  out_dir = "./results/trio_selection",
  date = format(Sys.Date(), "%Y%m%d"),
  threads = 8,
  dry_run = FALSE
)

usage <- function() {
  cat(
    "Usage: Rscript src/select/trio_selection.R --tree <file.nwk> [options]",
    "",
    "Options:",
    "  --tree FILE             Newick tree; leaf labels end in an NCBI accession. (required)",
    "  --idt-threshold N       Every pairwise identity must exceed this (default: 80).",
    "  --prescreen-margin N    Slack below the threshold when prescreening on tree",
    "                          distances, which only approximate aligned identity",
    "                          (default: 10).",
    "  --top-k N               Keep at most N outgroups per ingroup pair (default: 3).",
    "  --ingroup-min N         Require the ingroup pair's tree identity to be at least",
    "                          this before a trio becomes a candidate. On an ultrametric",
    "                          tree the outgroup is equidistant from both ingroups, so the",
    "                          ingroup pair's closeness is the only discriminating tree",
    "                          signal and hence the effective cost knob (default: 0 = off).",
    "  --min-contig-n50 BP     Drop a species whose best assembly has a contig N50 below",
    "                          this (Stage 0 quality gate). assembly_level and",
    "                          refseq_category mislabel fragmented assemblies, so contig",
    "                          N50 is the honest contiguity signal. Default 0 = off (keep",
    "                          every species); this only gates the on-tree assembly, so it",
    "                          is opt-in until the all-assembly redesign lands (e.g. pass",
    "                          1000000 for a primate tree).",
    "  --genome-dir DIR        Genome storage (default: ./genomes).",
    "  --out-dir DIR           Output and last-train cache (default: ./results/trio_selection).",
    "  --date YYYYMMDD         Run date, used in cache filenames (default: today).",
    "  --threads N             Threads for lastdb/last-train (default: 8).",
    "  --dry-run               Stop after Stage 1; download nothing and run no last-train.",
    "",
    sep = "\n"
  )
}

parse_args <- function(argv) {
  opts <- defaults
  i <- 1
  while (i <= length(argv)) {
    key <- argv[i]
    take <- function() {
      if (i + 1 > length(argv)) stop("Missing value for ", key, call. = FALSE)
      argv[i + 1]
    }
    switch(key,
      "--tree"             = { opts$tree <- take(); i <- i + 2 },
      "--idt-threshold"    = { opts$idt_threshold <- as.numeric(take()); i <- i + 2 },
      "--prescreen-margin" = { opts$prescreen_margin <- as.numeric(take()); i <- i + 2 },
      "--top-k"            = { opts$top_k <- as.integer(take()); i <- i + 2 },
      "--ingroup-min"      = { opts$ingroup_min <- as.numeric(take()); i <- i + 2 },
      "--min-contig-n50"   = { opts$min_contig_n50 <- as.numeric(take()); i <- i + 2 },
      "--genome-dir"       = { opts$genome_dir <- take(); i <- i + 2 },
      "--out-dir"          = { opts$out_dir <- take(); i <- i + 2 },
      "--date"             = { opts$date <- take(); i <- i + 2 },
      "--threads"          = { opts$threads <- as.integer(take()); i <- i + 2 },
      "--dry-run"          = { opts$dry_run <- TRUE; i <- i + 1 },
      "-h"                 = { usage(); quit(status = 0) },
      "--help"             = { usage(); quit(status = 0) },
      stop("Unknown option: ", key, call. = FALSE)
    )
  }
  if (is.null(opts$tree)) {
    usage()
    stop("--tree is required.", call. = FALSE)
  }
  opts
}

# --- helpers ----------------------------------------------------------------

log_stage <- function(...) cat("[trio_selection]", ..., "\n")

run <- function(command, args) {
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = ""))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop("Command failed (", status, "): ", command, " ", paste(args, collapse = " "),
         call. = FALSE)
  }
  output
}

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
             colClasses = "character")
}

# Leaf labels look like Genus_species_GCA_000000000.1; the accession is the tail.
ACCESSION_RE <- "GC[AF]_[0-9]+\\.[0-9]+$"

accession_from_label <- function(labels) {
  found <- regexpr(ACCESSION_RE, labels) > 0
  result <- rep(NA_character_, length(labels))
  result[found] <- regmatches(labels, regexpr(ACCESSION_RE, labels))
  result
}

# Short names must be stable per species and must NOT carry a trio slot number:
# the last-train cache is keyed on the species pair, so the same species has to
# produce the same name no matter which slot it lands in.  (sbst_fromDwl.sh's
# make_short_name appends the slot, which is why it cannot be reused here.)
make_short_names <- function(species) {
  base <- vapply(species, function(name) {
    tokens <- strsplit(name, "[ _]+")[[1]]
    first <- substr(tokens[1], 1, 3)
    second <- if (length(tokens) >= 2) substr(tokens[2], 1, 3) else substr(tokens[1], 4, 6)
    paste0(first, second)
  }, character(1), USE.NAMES = FALSE)

  # Disambiguate collisions so two species never share one cache entry.
  for (dup in unique(base[duplicated(base)])) {
    idx <- which(base == dup)
    base[idx] <- paste0(dup, seq_along(idx))
  }
  base
}

descendant_tips <- function(tree, node) {
  kids <- tree$edge[tree$edge[, 1] == node, 2]
  tips <- kids[kids <= Ntip(tree)]
  internal <- kids[kids > Ntip(tree)]
  out <- tips
  for (u in internal) out <- c(out, descendant_tips(tree, u))
  out
}

sister_tips <- function(tree, node) {
  parent <- tree$edge[tree$edge[, 2] == node, 1]
  if (!length(parent)) return(integer(0))
  siblings <- setdiff(tree$edge[tree$edge[, 1] == parent, 2], node)
  unlist(lapply(siblings, descendant_tips, tree = tree))
}

# --- Stage 0: prune the tree to one assembly per species ---------------------

prune_to_best_assembly <- function(tree, out_dir, min_contig_n50 = 0) {
  accessions <- accession_from_label(tree$tip.label)
  if (anyNA(accessions)) {
    stop("Leaf label(s) without a trailing NCBI accession: ",
         paste(head(tree$tip.label[is.na(accessions)], 5), collapse = ", "), call. = FALSE)
  }

  acc_file <- file.path(out_dir, "tree_accessions.txt")
  meta_file <- file.path(out_dir, "assembly_metadata.tsv")
  writeLines(unique(accessions), acc_file)

  log_stage("Stage 0: fetching NCBI metadata for", length(unique(accessions)), "accessions")
  run("python3", c(shQuote(FETCH_METADATA),
                   "--accessions", shQuote(acc_file),
                   "--out", shQuote(meta_file)))

  meta <- read_tsv(meta_file)
  meta$contig_n50 <- suppressWarnings(as.numeric(meta$contig_n50))

  leaves <- tibble::tibble(tip = tree$tip.label, accession = accessions) %>%
    inner_join(meta, by = "accession")

  dropped <- setdiff(tree$tip.label, leaves$tip)
  if (length(dropped)) {
    log_stage("  dropping", length(dropped), "leaf/leaves with no NCBI record:",
              paste(dropped, collapse = ", "))
  }

  best <- leaves %>%
    mutate(
      is_reference = as.integer(refseq_category == "reference genome"),
      level_rank = ifelse(
        assembly_level %in% names(ASSEMBLY_LEVEL_RANK),
        ASSEMBLY_LEVEL_RANK[assembly_level],
        length(ASSEMBLY_LEVEL_RANK) + 1L
      ),
      contig_n50 = ifelse(is.na(contig_n50), 0, contig_n50)
    ) %>%
    arrange(species, desc(is_reference), level_rank, desc(contig_n50)) %>%
    group_by(species) %>%
    slice_head(n = 1) %>%
    ungroup()

  log_stage(" ", nrow(leaves), "assemblies ->", nrow(best), "species")

  # Species-level quality gate.  assembly_level labels every primate assembly
  # "Chromosome" even when its contig N50 is ~13 kb, and refseq_category
  # "reference genome" likewise covers badly fragmented assemblies, so neither is a
  # usable quality signal -- gate on contig N50, the honest contiguity measure.  A
  # species whose best assembly is below the floor is dropped entirely (all of its
  # trios go with it), because a fragmented genome shrinks the alignable and
  # callable-CDS fraction and biases which genes and sequence contexts are observed
  # (contig N50 measures continuity only -- not completeness, base accuracy, or
  # contamination, and can even be inflated by mis-joins).
  #
  # LIMITATION: this gates only the assembly that happens to be on the tree, so a
  # species whose tree accession is fragmented is dropped even when a better current
  # assembly exists off-tree (e.g. Pongo abelii: tree 16 kb vs current best 146 Mb).
  # The planned fix is to fetch every current assembly per species, gate the set, and
  # drop the species only when no candidate passes.  contig N50 is also assembly-size
  # relative; NG50 (against expected genome size) is preferable across distant taxa.
  if (min_contig_n50 > 0) {
    low <- best[best$contig_n50 < min_contig_n50, , drop = FALSE]
    if (nrow(low)) {
      log_stage("  quality gate: dropping", nrow(low), "species with contig N50 <",
                format(min_contig_n50, scientific = FALSE, big.mark = ","), "bp:",
                paste(sprintf("%s (%.3g)", low$species, low$contig_n50), collapse = ", "))
      best <- best[best$contig_n50 >= min_contig_n50, , drop = FALSE]
    }
    log_stage("  ", nrow(best), "species pass the contig-N50 gate")
  }

  pruned <- drop.tip(tree, setdiff(tree$tip.label, best$tip))
  best <- best[match(pruned$tip.label, best$tip), ]
  best$short_name <- make_short_names(best$species)

  list(tree = pruned, leaves = best)
}

# --- Stage 1: enumerate trios and prescreen on the tree ----------------------

enumerate_trios <- function(tree, leaves) {
  # Branch lengths are percent sequence difference, so identity is 100 - distance.
  tree_identity <- 100 - cophenetic(tree)
  tips <- tree$tip.label
  info <- leaves[match(tips, leaves$tip), ]

  rows <- list()
  for (node in (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)) {
    ingroup <- descendant_tips(tree, node)
    if (length(ingroup) < 2) next
    outgroup <- sister_tips(tree, node)
    if (!length(outgroup)) next

    pairs <- combn(ingroup, 2)
    for (p in seq_len(ncol(pairs))) {
      in1 <- pairs[1, p]
      in2 <- pairs[2, p]
      for (out in outgroup) {
        rows[[length(rows) + 1]] <- tibble::tibble(
          mrca_node = node,
          out_tip = tips[out], in1_tip = tips[in1], in2_tip = tips[in2],
          out_acc = info$accession[out],
          in1_acc = info$accession[in1],
          in2_acc = info$accession[in2],
          out_species = info$species[out],
          in1_species = info$species[in1],
          in2_species = info$species[in2],
          out_short = info$short_name[out],
          in1_short = info$short_name[in1],
          in2_short = info$short_name[in2],
          # Slot 1 is the outgroup, matching org1/org2/org3 in the rest of the pipeline.
          genus_1 = info$genus[out],
          genus_2 = info$genus[in1],
          genus_3 = info$genus[in2],
          tree_idt_12 = tree_identity[tips[out], tips[in1]],
          tree_idt_13 = tree_identity[tips[out], tips[in2]],
          tree_idt_23 = tree_identity[tips[in1], tips[in2]]
        )
      }
    }
  }
  bind_rows(rows)
}

prescreen_trios <- function(trios, opts) {
  # A trio whose genera split two-vs-one -- the outgroup congeneric with exactly
  # one ingroup species -- is excluded by the thesis rule whatever the identities
  # say, so it can be discarded before spending a single download.
  floor_idt <- opts$idt_threshold - opts$prescreen_margin

  trios %>%
    mutate(
      genus_two_vs_one =
        (genus_1 == genus_2 & genus_1 != genus_3) |
        (genus_1 == genus_3 & genus_1 != genus_2),
      tree_idt_min = pmin(tree_idt_12, tree_idt_13, tree_idt_23),
      outgroup_idt = pmin(tree_idt_12, tree_idt_13)
    ) %>%
    filter(!genus_two_vs_one, tree_idt_min >= floor_idt, tree_idt_23 >= opts$ingroup_min) %>%
    # The closest admissible outgroup resolves substitutions best, so rank on it.
    group_by(in1_tip, in2_tip) %>%
    arrange(desc(outgroup_idt), .by_group = TRUE) %>%
    slice_head(n = opts$top_k) %>%
    ungroup()
}

# --- Stage 2: download the surviving accessions ------------------------------

download_genomes <- function(trios, opts) {
  wanted <- unique(c(trios$out_acc, trios$in1_acc, trios$in2_acc))
  log_stage("Stage 2: downloading", length(wanted), "genomes into", opts$genome_dir)

  fasta <- character(0)
  for (accession in wanted) {
    result <- run("bash", c(shQuote(DWL_ORGANISM), accession,
                            "--out-dir", shQuote(opts$genome_dir)))
    # dwl_organism.sh prints dir_name|fasta_path|summary_json|raw_name|ncbi_full|tax_json
    fields <- strsplit(tail(result, 1), "|", fixed = TRUE)[[1]]
    if (length(fields) < 2 || !nzchar(fields[2])) {
      stop("dwl_organism.sh returned no FASTA path for ", accession, call. = FALSE)
    }
    # last_train.sh cds into its lastdb directory, so the FASTA path must be absolute.
    fasta[accession] <- normalizePath(fields[2], mustWork = TRUE)
  }
  fasta
}

# --- Stage 3: last-train once per unique species pair ------------------------

unique_species_pairs <- function(trios) {
  pairs <- bind_rows(
    trios %>% transmute(a_acc = out_acc, a_short = out_short,
                        b_acc = in1_acc, b_short = in1_short),
    trios %>% transmute(a_acc = out_acc, a_short = out_short,
                        b_acc = in2_acc, b_short = in2_short),
    trios %>% transmute(a_acc = in1_acc, a_short = in1_short,
                        b_acc = in2_acc, b_short = in2_short)
  )
  # Canonical (unordered) orientation: one species pair must map to one cache
  # entry regardless of which slots it occupied in any given trio.
  pairs %>%
    mutate(
      first_acc = ifelse(a_short <= b_short, a_acc, b_acc),
      first_short = ifelse(a_short <= b_short, a_short, b_short),
      second_acc = ifelse(a_short <= b_short, b_acc, a_acc),
      second_short = ifelse(a_short <= b_short, b_short, a_short)
    ) %>%
    distinct(first_acc, first_short, second_acc, second_short)
}

train_file_path <- function(cache_dir, short_a, short_b, date) {
  file.path(cache_dir, sprintf("%s2%s_%s.train", short_a, short_b, date))
}

parse_train_identity <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  line <- grep(IDENTITY_PREFIX, readLines(path, warn = FALSE), fixed = TRUE, value = TRUE)
  if (!length(line)) return(NA_real_)
  value <- regmatches(line[1], regexpr("[0-9]+\\.?[0-9]*", line[1]))
  if (!length(value)) NA_real_ else as.numeric(value)
}

train_species_pairs <- function(pairs, fasta, opts) {
  cache_dir <- file.path(opts$out_dir, "train_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(THREAD_NUM = opts$threads)

  log_stage("Stage 3: last-train on", nrow(pairs), "unique species pairs")

  identity <- numeric(nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    short_a <- pairs$first_short[i]
    short_b <- pairs$second_short[i]
    train <- train_file_path(cache_dir, short_a, short_b, opts$date)

    if (file.exists(train)) {
      log_stage("  cached:", basename(train))
    } else {
      log_stage("  training:", short_a, "vs", short_b)
      # last_train.sh is idempotent: it skips lastdb and last-train when their
      # outputs already exist, so a shared cache_dir gives cross-trio reuse.
      run("bash", c(shQuote(LAST_TRAIN), opts$date,
                    shQuote(fasta[[pairs$first_acc[i]]]),
                    shQuote(fasta[[pairs$second_acc[i]]]),
                    short_a, short_b, shQuote(cache_dir)))
    }
    identity[i] <- parse_train_identity(train)
  }

  pairs$identity <- identity
  pairs
}

attach_identities <- function(trios, pairs) {
  lookup <- setNames(pairs$identity, paste(pairs$first_short, pairs$second_short, sep = "|"))
  key <- function(x, y) ifelse(x <= y, paste(x, y, sep = "|"), paste(y, x, sep = "|"))

  trios %>%
    mutate(
      idt_12 = unname(lookup[key(out_short, in1_short)]),
      idt_13 = unname(lookup[key(out_short, in2_short)]),
      idt_23 = unname(lookup[key(in1_short, in2_short)])
    )
}

# --- Stage 4: apply the thesis rule ------------------------------------------

apply_thesis_rule <- function(trios, opts) {
  trios_file <- file.path(opts$out_dir, "candidate_trios.tsv")
  verdict_file <- file.path(opts$out_dir, "trio_verdicts.tsv")

  write_tsv(trios, trios_file)
  log_stage("Stage 4: applying the thesis rule via trio_filter.py")
  run("python3", c(shQuote(TRIO_FILTER),
                   "--trios", shQuote(trios_file),
                   "--idt-threshold", opts$idt_threshold,
                   "--out", shQuote(verdict_file)))

  read_tsv(verdict_file)
}

# --- main --------------------------------------------------------------------

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

  tree <- read.tree(opts$tree)
  log_stage("tree:", opts$tree, "-", Ntip(tree), "leaves")

  pruned <- prune_to_best_assembly(tree, opts$out_dir, opts$min_contig_n50)
  write_tsv(pruned$leaves, file.path(opts$out_dir, "selected_assemblies.tsv"))

  trios <- enumerate_trios(pruned$tree, pruned$leaves)
  log_stage("Stage 1:", nrow(trios), "trios enumerated")

  candidates <- prescreen_trios(trios, opts)
  log_stage(" ", nrow(candidates), "candidates after prescreen",
            sprintf("(two-vs-one genera dropped; tree identity >= %g; ingroup identity >= %g; top %d outgroups per ingroup pair)",
                    opts$idt_threshold - opts$prescreen_margin, opts$ingroup_min, opts$top_k))

  if (!nrow(candidates)) {
    log_stage("No candidate trios survived the prescreen; nothing to do.")
    return(invisible(NULL))
  }

  if (opts$dry_run) {
    out <- file.path(opts$out_dir, "candidate_trios.tsv")
    write_tsv(candidates, out)
    log_stage("--dry-run: stopping after Stage 1. Candidates written to", out)
    return(invisible(NULL))
  }

  fasta <- download_genomes(candidates, opts)
  pairs <- train_species_pairs(unique_species_pairs(candidates), fasta, opts)
  candidates <- attach_identities(candidates, pairs)

  verdicts <- apply_thesis_rule(candidates, opts)
  selected <- verdicts[verdicts$excluded == "FALSE", , drop = FALSE]

  selected_file <- file.path(opts$out_dir, "selected_trios.tsv")
  write_tsv(selected, selected_file)

  log_stage("Selected", nrow(selected), "of", nrow(verdicts), "trios ->", selected_file)
}

if (sys.nframe() == 0) {
  main()
}
