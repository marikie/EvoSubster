#!/usr/bin/env Rscript

# Select substitution-analysis trios from a Newick tree.
#
# Original script written by Ryota Ishii; rewritten to apply the thesis
# (section 2.3) rule -- measured sequence percent identity plus genus names --
# instead of tree distances.
#
# One outgroup per ingroup pair.  The tree is used only for TOPOLOGY: for an
# ingroup couple it gives the candidate outgroups ordered nearest -> farthest
# (sister of their MRCA, then sister of the grandparent, ...).  The tree is
# ultrametric, so its distances cannot rank outgroups -- the real decision needs
# last-train identity.  For each pair we walk the candidates nearest-first,
# download + last-train on demand, and STOP at the first outgroup that passes the
# thesis rule.  The strict thesis ordering (idt_12 < idt_23 AND idt_13 < idt_23)
# is itself the "not too close" boundary, and --idt-threshold is the "not too far"
# boundary, so nearest-first + first-pass yields the closest still-external
# outgroup while training as few pairs as possible.
#
# Stage 0  prune the tree to one assembly per species      [NCBI metadata, no download]
# Stage 1  per ingroup pair: search outgroups near->far,   [dwl_organism.sh + last_train.sh
#          stop at the first thesis-rule pass -> one trio    on demand; trio_filter.py to judge]
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
  max_outgroup_tries = 5,
  ingroup_pairing = "matching",
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
    "  --max-outgroup-tries N  Give up on an ingroup pair after training this many outgroup",
    "                          candidates without a thesis-rule pass (default: 5). Caps the",
    "                          per-pair last-train cost.",
    "  --ingroup-pairing MODE  Which ingroup couples to consider (default: matching).",
    "                            matching  greedy closest-first disjoint matching -- each",
    "                                      species used at most once, ~n/2 independent",
    "                                      sister-pair couples (phylogenetic-pairs design).",
    "                            all       every tip pair (exhaustive; C(n,2) couples).",
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
    "  --dry-run               Enumerate ingroup pairs + candidate-outgroup order only;",
    "                          download nothing and run no last-train.",
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
      "--tree"               = { opts$tree <- take(); i <- i + 2 },
      "--idt-threshold"      = { opts$idt_threshold <- as.numeric(take()); i <- i + 2 },
      "--max-outgroup-tries" = { opts$max_outgroup_tries <- as.integer(take()); i <- i + 2 },
      "--ingroup-pairing"    = { opts$ingroup_pairing <- take(); i <- i + 2 },
      "--min-contig-n50"     = { opts$min_contig_n50 <- as.numeric(take()); i <- i + 2 },
      "--genome-dir"         = { opts$genome_dir <- take(); i <- i + 2 },
      "--out-dir"            = { opts$out_dir <- take(); i <- i + 2 },
      "--date"               = { opts$date <- take(); i <- i + 2 },
      "--threads"            = { opts$threads <- as.integer(take()); i <- i + 2 },
      "--dry-run"            = { opts$dry_run <- TRUE; i <- i + 1 },
      "-h"                   = { usage(); quit(status = 0) },
      "--help"               = { usage(); quit(status = 0) },
      stop("Unknown option: ", key, call. = FALSE)
    )
  }
  if (is.null(opts$tree)) {
    usage()
    stop("--tree is required.", call. = FALSE)
  }
  if (!opts$ingroup_pairing %in% c("matching", "all")) {
    stop("--ingroup-pairing must be 'matching' or 'all' (got: ", opts$ingroup_pairing, ")",
         call. = FALSE)
  }
  opts
}

# --- helpers ----------------------------------------------------------------

log_stage <- function(...) cat("[trio_selection]", ..., "\n")

run <- function(command, args, quiet_stderr = FALSE) {
  output <- suppressWarnings(system2(command, args, stdout = TRUE,
                                     stderr = if (quiet_stderr) FALSE else ""))
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
  # A tip is its own only descendant tip; otherwise recurse into the children.
  # (Without the tip base case a single-leaf sister clade is silently dropped, so
  # a pair whose only outgroups are lone leaves would get no candidates at all.)
  if (node <= Ntip(tree)) return(node)
  kids <- tree$edge[tree$edge[, 1] == node, 2]
  unlist(lapply(kids, descendant_tips, tree = tree))
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

# --- Stage 1: per ingroup pair, ordered candidate outgroups (near -> far) -----

# For the ingroup couple (tipA, tipB), list candidate outgroups nearest first:
# the sister clade of their MRCA, then the sister of the grandparent, ... up to
# the root.  On an ultrametric tree the members of one sister clade are tied, so
# they are ordered by Stage-0 assembly quality (contig N50) as a tie-break.
# Returns an empty vector when the pair's MRCA is the root (no outgroup exists).
candidate_outgroups <- function(tree, tipA, tipB, leaves) {
  mrca <- getMRCA(tree, c(tipA, tipB))
  root <- Ntip(tree) + 1L
  if (is.null(mrca)) return(character(0))

  cands <- integer(0)
  node <- mrca
  while (length(node) == 1 && node != root) {
    sibs <- sister_tips(tree, node)
    if (length(sibs)) {
      q <- leaves$contig_n50[match(tree$tip.label[sibs], leaves$tip)]
      q[is.na(q)] <- 0
      cands <- c(cands, sibs[order(-q)])
    }
    parent <- tree$edge[tree$edge[, 2] == node, 1]
    node <- if (length(parent)) parent else root
  }
  tree$tip.label[cands]
}

# Choose which ingroup couples to consider, using tree topology only (no last-train).
#
#   "matching" (default): a phylogenetically-independent sister-pair design -- pair each
#     species with its closest relative and use each species at most once.  Candidate
#     couples are ranked closest-first by TOPOLOGICAL distance (the number of edges between
#     the tips: every edge counts as 1, so a cherry = distance 2 and beats any deeper pair)
#     and greedily matched, giving ~n/2 disjoint, independent sister-pair couples.  Topology
#     (not branch length) is used deliberately: rate variation or additive (substitutions/
#     site) branch lengths could otherwise rank a non-sister pair with short branches ahead
#     of a true cherry with long branches.  Only couples with an available outgroup (MRCA is
#     not the root) are eligible, so a species is never spent on a couple that could yield no
#     trio.  A species left unmatched on a ladder-like tree simply produces no trio.
#
#   "all": every tip pair (exhaustive C(n,2)); the pre-restriction behavior, kept as an
#     opt-out / validation path.
#
# Returns a list of c(tipA, tipB) character pairs.
ingroup_pairs <- function(tree, leaves, mode = "matching") {
  tips <- tree$tip.label
  n <- length(tips)
  if (n < 2) return(list())

  if (mode == "all") {
    combs <- utils::combn(n, 2)
    return(lapply(seq_len(ncol(combs)), function(k) tips[combs[, k]]))
  }

  # Closeness is measured on the tree TOPOLOGY, not branch lengths: set every edge to 1
  # so a pair's distance is the number of edges between them.  A cherry (true sisters) is
  # always distance 2 -- the minimum -- so cherries are matched before any non-sister pair.
  # Real branch lengths would break this: on a tree with rate variation or additive
  # (substitutions/site) lengths a non-sister pair with short branches can have a smaller
  # patristic distance than a true sister pair with long branches, so branch-length ranking
  # could pair non-sisters.  Topology also makes the matching independent of the tree's units
  # (time vs substitutions) and of whether it is ultrametric.
  t2 <- tree
  t2$edge.length <- rep(1, nrow(t2$edge))
  D <- cophenetic(t2)[tips, tips]
  M <- mrca(tree)
  root <- Ntip(tree) + 1L

  ij <- which(upper.tri(D), arr.ind = TRUE)     # candidate i<j pairs (by tip index)
  eligible <- M[ij] != root                     # keep only couples that have an outgroup
  ij <- ij[eligible, , drop = FALSE]
  if (!nrow(ij)) return(list())

  # Deterministic ordering: fewest edges (closest / most-nested) first, then better combined
  # assembly quality, then tip names -- so the same tree always yields the same matching.
  q <- leaves$contig_n50[match(tips, leaves$tip)]
  q[is.na(q)] <- 0
  score <- q[ij[, 1]] + q[ij[, 2]]
  ord <- order(D[ij], -score, tips[ij[, 1]], tips[ij[, 2]])

  used <- logical(n)
  pairs <- list()
  for (k in ord) {
    a <- ij[k, 1]; b <- ij[k, 2]
    if (used[a] || used[b]) next
    used[a] <- TRUE; used[b] <- TRUE
    pairs[[length(pairs) + 1L]] <- c(tips[a], tips[b])
  }
  pairs
}

# Enumerate the selected ingroup couples with their ordered candidate-outgroup list.
# Used by --dry-run to expose the search order without training.
enumerate_pairs_dry <- function(tree, leaves, mode = "matching") {
  rows <- list()
  for (p in ingroup_pairs(tree, leaves, mode)) {
    cands <- candidate_outgroups(tree, p[1], p[2], leaves)
    if (!length(cands)) next
    rows[[length(rows) + 1]] <- data.frame(
      in1_tip = p[1], in2_tip = p[2],
      n_candidates = length(cands),
      candidate_order = paste(cands, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Mirror of trio_filter.py pattern5: the outgroup is congeneric with exactly one
# ingroup species.  The thesis rule excludes such a trio whatever the identities
# say, so skipping it here avoids a wasted last-train.
is_two_vs_one <- function(g_out, g_in1, g_in2) {
  o <- tolower(g_out); a <- tolower(g_in1); b <- tolower(g_in2)
  xor(o == a, o == b)
}

# --- on-demand download + last-train, memoized ------------------------------

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

# Build closures that download a genome and last-train a species pair at most
# once each (memoized), so a species trained for one ingroup pair is reused for
# every other.  The shared cache dir also lets last_train.sh skip work across runs.
make_fetchers <- function(leaves, opts) {
  fasta_cache <- new.env(parent = emptyenv())
  idt_cache <- new.env(parent = emptyenv())
  counters <- new.env(parent = emptyenv())
  counters$downloads <- 0L
  counters$trains <- 0L

  cache_dir <- file.path(opts$out_dir, "train_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(THREAD_NUM = opts$threads)

  get_fasta <- function(acc) {
    hit <- fasta_cache[[acc]]
    if (!is.null(hit)) return(hit)
    result <- run("bash", c(shQuote(DWL_ORGANISM), acc, "--out-dir", shQuote(opts$genome_dir)))
    # dwl_organism.sh prints dir_name|fasta_path|summary_json|raw_name|ncbi_full|tax_json
    fields <- strsplit(tail(result, 1), "|", fixed = TRUE)[[1]]
    if (length(fields) < 2 || !nzchar(fields[2])) {
      stop("dwl_organism.sh returned no FASTA path for ", acc, call. = FALSE)
    }
    # last_train.sh cds into its lastdb directory, so the FASTA path must be absolute.
    path <- normalizePath(fields[2], mustWork = TRUE)
    fasta_cache[[acc]] <- path
    counters$downloads <- counters$downloads + 1L
    path
  }

  get_identity <- function(short_a, acc_a, short_b, acc_b) {
    # Canonical (unordered) orientation so one species pair is one cache entry.
    if (short_a <= short_b) {
      s1 <- short_a; a1 <- acc_a; s2 <- short_b; a2 <- acc_b
    } else {
      s1 <- short_b; a1 <- acc_b; s2 <- short_a; a2 <- acc_a
    }
    key <- paste(s1, s2, sep = "|")
    hit <- idt_cache[[key]]
    if (!is.null(hit)) return(hit)

    train <- train_file_path(cache_dir, s1, s2, opts$date)
    if (!file.exists(train)) {
      fa1 <- get_fasta(a1); fa2 <- get_fasta(a2)
      run("bash", c(shQuote(LAST_TRAIN), opts$date, shQuote(fa1), shQuote(fa2),
                    s1, s2, shQuote(cache_dir)))
      counters$trains <- counters$trains + 1L
    }
    val <- parse_train_identity(train)
    idt_cache[[key]] <- val
    val
  }

  list(get_identity = get_identity, counters = counters)
}

# --- the thesis rule, judged one candidate at a time ------------------------

# Shell one candidate trio through trio_filter.py (the single source of the
# thesis rule) and return its verdict row.
evaluate_one <- function(row, opts, tmp_in, tmp_out) {
  write_tsv(row, tmp_in)
  run("python3", c(shQuote(TRIO_FILTER),
                   "--trios", shQuote(tmp_in),
                   "--idt-threshold", opts$idt_threshold,
                   "--out", shQuote(tmp_out)),
      quiet_stderr = TRUE)
  read_tsv(tmp_out)[1, , drop = FALSE]
}

build_row <- function(infoO, infoA, infoB, idt_12, idt_13, idt_23, tries) {
  data.frame(
    out_tip = infoO$tip, in1_tip = infoA$tip, in2_tip = infoB$tip,
    out_acc = infoO$accession, in1_acc = infoA$accession, in2_acc = infoB$accession,
    out_species = infoO$species, in1_species = infoA$species, in2_species = infoB$species,
    out_short = infoO$short_name, in1_short = infoA$short_name, in2_short = infoB$short_name,
    # Slot 1 is the outgroup, matching org1/org2/org3 in the rest of the pipeline.
    genus_1 = infoO$genus, genus_2 = infoA$genus, genus_3 = infoB$genus,
    idt_12 = idt_12, idt_13 = idt_13, idt_23 = idt_23,
    outgroup_tries = tries,
    stringsAsFactors = FALSE
  )
}

# --- Stage 1 driver: one nearest-but-external outgroup per ingroup pair -------

select_trios <- function(tree, leaves, opts, fetch = NULL) {
  if (is.null(fetch)) fetch <- make_fetchers(leaves, opts)
  get_identity <- fetch$get_identity
  threshold <- opts$idt_threshold
  tmp_in <- file.path(opts$out_dir, ".candidate_in.tsv")
  tmp_out <- file.path(opts$out_dir, ".candidate_out.tsv")

  info_of <- function(tip) leaves[match(tip, leaves$tip), , drop = FALSE]

  couples <- ingroup_pairs(tree, leaves, opts$ingroup_pairing)
  covered <- unique(unlist(couples))
  uncovered <- setdiff(tree$tip.label, covered)
  results <- list()
  no_outgroup <- 0L
  too_diverged <- 0L
  no_valid_outgroup <- 0L

  log_stage(sprintf("Stage 1: %s pairing -> %d ingroup couple(s); %d species not paired",
                    opts$ingroup_pairing, length(couples), length(uncovered)))
  log_stage("Stage 1: searching one outgroup per ingroup pair (nearest-first, first-pass)")

  for (couple in couples) {
    tipA <- couple[1]; tipB <- couple[2]
    cands <- candidate_outgroups(tree, tipA, tipB, leaves)
    if (!length(cands)) { no_outgroup <- no_outgroup + 1L; next }

    infoA <- info_of(tipA); infoB <- info_of(tipB)
    idt_23 <- get_identity(infoA$short_name, infoA$accession,
                           infoB$short_name, infoB$accession)
    if (is.na(idt_23) || idt_23 <= threshold) { too_diverged <- too_diverged + 1L; next }

    chosen <- NULL
    tries <- 0L
    for (og in cands) {
      infoO <- info_of(og)
      # Free pre-filter: a two-vs-one genus split is always excluded -> never train it.
      if (is_two_vs_one(infoO$genus, infoA$genus, infoB$genus)) next
      if (tries >= opts$max_outgroup_tries) break
      tries <- tries + 1L

      idt_12 <- get_identity(infoO$short_name, infoO$accession,
                             infoA$short_name, infoA$accession)
      idt_13 <- get_identity(infoO$short_name, infoO$accession,
                             infoB$short_name, infoB$accession)
      if (is.na(idt_12) || is.na(idt_13)) next  # unusable alignment, try the next candidate

      row <- build_row(infoO, infoA, infoB, idt_12, idt_13, idt_23, tries)
      verdict <- evaluate_one(row, opts, tmp_in, tmp_out)

      if (identical(verdict$excluded, "FALSE")) {           # pass -> fix this trio, stop
        chosen <- verdict
        break
      }
      if (identical(verdict$idt_threshold_condition, "FALSE")) {
        # too far: an outgroup-ingroup identity fell below the threshold, and every
        # farther candidate is more diverged still -> the window has closed.
        break
      }
      # else "too close" (strict ordering failed) -> the next candidate is farther out.
    }

    if (!is.null(chosen)) {
      results[[length(results) + 1]] <- chosen
      log_stage(sprintf("  %s + %s  ->  %s  (idt 12/13/23 = %s/%s/%s; %d tried)",
                        infoA$species, infoB$species, chosen$out_species,
                        chosen$idt_12, chosen$idt_13, chosen$idt_23, tries))
    } else {
      no_valid_outgroup <- no_valid_outgroup + 1L
    }
  }

  log_stage(sprintf("last-train runs: %d | downloads: %d", fetch$counters$trains,
                    fetch$counters$downloads))
  log_stage(sprintf("ingroup couples: %d selected | %d no outgroup on tree | %d couple below threshold | %d no external outgroup found | %d species unpaired",
                    length(results), no_outgroup, too_diverged, no_valid_outgroup,
                    length(uncovered)))

  if (!length(results)) return(data.frame())
  do.call(rbind, results)
}

# --- main --------------------------------------------------------------------

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

  tree <- read.tree(opts$tree)
  log_stage("tree:", opts$tree, "-", Ntip(tree), "leaves")

  pruned <- prune_to_best_assembly(tree, opts$out_dir, opts$min_contig_n50)
  write_tsv(pruned$leaves, file.path(opts$out_dir, "selected_assemblies.tsv"))
  tree <- pruned$tree
  leaves <- pruned$leaves

  if (opts$dry_run) {
    pairs <- enumerate_pairs_dry(tree, leaves, opts$ingroup_pairing)
    out <- file.path(opts$out_dir, "candidate_outgroups.tsv")
    write_tsv(pairs, out)
    log_stage("Stage 1 (dry-run):", nrow(pairs), "ingroup couple(s) with an available",
              "outgroup (", opts$ingroup_pairing, "pairing); candidate order ->", out)
    return(invisible(NULL))
  }

  selected <- select_trios(tree, leaves, opts)
  selected_file <- file.path(opts$out_dir, "selected_trios.tsv")
  write_tsv(selected, selected_file)
  log_stage("Selected", nrow(selected), "trios (one per ingroup pair) ->", selected_file)
}

if (sys.nframe() == 0) {
  main()
}
