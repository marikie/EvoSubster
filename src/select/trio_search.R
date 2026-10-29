# Helpers loaded by trio_selection.R into its calling environment.

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
      # Relative (size-normalized) contig N50 tie-break: compare contiguity fairly across
      # species with very different genome sizes, not raw base counts.
      q <- leaves$rel_contig_n50[match(tree$tip.label[sibs], leaves$tip)]
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
  # assembly quality (relative contig N50, size-normalized), then tip names -- so the same
  # tree always yields the same matching.
  q <- leaves$rel_contig_n50[match(tips, leaves$tip)]
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

# --- the thesis rule, judged one candidate at a time ------------------------

# Shell one candidate trio through trio_filter.py (the single source of the
# thesis rule) and return its verdict row.
evaluate_one <- function(row, opts, tmp_in, tmp_out) {
  opts <- normalize_identity_options(opts)
  write_tsv(row, tmp_in)
  run("python3", c(shQuote(TRIO_FILTER),
                   "--trios", shQuote(tmp_in),
                   "--idt-threshold", opts$idt_threshold,
                   "--ingroup-idt-threshold", opts$ingroup_idt_threshold,
                   "--outgroup-idt-threshold", opts$outgroup_idt_threshold,
                   "--min-idt-gap-between-out-in", opts$min_idt_gap_between_out_in,
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

select_trio_for_pair <- function(tree, leaves, couple, opts, get_identity, tmp_in, tmp_out) {
  info_of <- function(tip) leaves[match(tip, leaves$tip), , drop = FALSE]
  tipA <- couple[1]; tipB <- couple[2]
  cands <- candidate_outgroups(tree, tipA, tipB, leaves)
  if (!length(cands)) return(list(status = "no_outgroup"))

  infoA <- info_of(tipA); infoB <- info_of(tipB)
  idt_23 <- get_identity(infoA$accession, infoB$accession)
  if (is.na(idt_23) || idt_23 <= opts$ingroup_idt_threshold) return(list(status = "too_diverged"))

  chosen <- NULL
  tries <- 0L
  for (og in cands) {
    infoO <- info_of(og)
    # Free pre-filter: a two-vs-one genus split is always excluded -> never train it.
    if (is_two_vs_one(infoO$genus, infoA$genus, infoB$genus)) next
    if (tries >= opts$max_outgroup_tries) break
    tries <- tries + 1L

    idt_12 <- get_identity(infoO$accession, infoA$accession)
    idt_13 <- get_identity(infoO$accession, infoB$accession)
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
    # else "too close" (configured gap failed) -> try the next, farther candidate.
  }
  list(status = if (is.null(chosen)) "no_valid_outgroup" else "selected",
       chosen = chosen, infoA = infoA, infoB = infoB, tries = tries)
}

select_trios <- function(tree, leaves, opts, fetch = NULL) {
  opts <- normalize_identity_options(opts)
  if (is.null(fetch)) fetch <- make_fetchers(leaves, opts)
  get_identity <- fetch$get_identity
  tmp_in <- file.path(opts$out_dir, ".candidate_in.tsv")
  tmp_out <- file.path(opts$out_dir, ".candidate_out.tsv")

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
    result <- select_trio_for_pair(tree, leaves, couple, opts, get_identity, tmp_in, tmp_out)
    if (result$status == "no_outgroup") { no_outgroup <- no_outgroup + 1L; next }
    if (result$status == "too_diverged") { too_diverged <- too_diverged + 1L; next }
    chosen <- result$chosen
    infoA <- result$infoA; infoB <- result$infoB; tries <- result$tries
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
