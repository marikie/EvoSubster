library(ape)
library(dplyr)
library(purrr)
library(tibble)

# --- 0) Load tree and compute patristic identities ---------------------------
tr <- read.tree("Eupercaria_exported.nwk")
patr_dist <- cophenetic(tr)
patr_idt  <- 100 - patr_dist   # patristic % identity

# thresholds
INGROUP_ID_MIN <- 95
INGROUP_ID_MAX <- 99
OUTGRP_ID_MIN  <- 85
OUTGRP_ID_MAX  <- 99
BALANCE_DELTA  <- 2

# --- 1) helpers --------------------------------------------------------------
descendant_tips <- function(tree, node) {
  kids <- tree$edge[tree$edge[,1] == node, 2]
  tips <- kids[kids <= Ntip(tree)]
  ints <- kids[kids >  Ntip(tree)]
  out  <- if (length(tips)) tips else integer(0)
  for (u in ints) out <- c(out, descendant_tips(tree, u))
  out
}

sister_tips <- function(tree, node) {
  parent <- tree$edge[tree$edge[,2] == node, 1]
  if (!length(parent)) return(integer(0))
  sibs <- tree$edge[tree$edge[,1] == parent, 2]
  sibs <- setdiff(sibs, node)
  unlist(lapply(sibs, descendant_tips, tree = tree))
}

parsimony_ratio <- function(idt_AB, idt_AC, idt_BC) {
  pA <- (100 - idt_AB) / 2
  pB <- pA
  idt_out <- (idt_AC + idt_BC) / 2
  pC <- 100 - idt_out - pA
  ratio <- ((100 - pC) / pC) * 3
  bad <- !is.finite(ratio) | is.na(pC) | pC <= 0 | pC >= 100
  ratio[bad] <- NA_real_
  ratio
}

# --- 2) enumerate trios -----------------------------------------------------
tip_labels <- tr$tip.label
internal_nodes <- (Ntip(tr) + 1):(Ntip(tr) + tr$Nnode)

trios <- list()

for (u in internal_nodes) {
  ing_tips <- descendant_tips(tr, u)
  if (length(ing_tips) < 2) next
  out_tips <- sister_tips(tr, u)
  if (!length(out_tips)) next
  
  ing_pairs <- t(combn(ing_tips, 2))
  for (r in seq_len(nrow(ing_pairs))) {
    i <- ing_pairs[r,1]; j <- ing_pairs[r,2]
    for (k in out_tips) {
      A <- tip_labels[i]; B <- tip_labels[j]; C <- tip_labels[k]
      idt_AB <- patr_idt[A, B]
      idt_AC <- patr_idt[A, C]
      idt_BC <- patr_idt[B, C]
      
      trios[[length(trios)+1]] <- tibble(
        ing1 = A, ing2 = B, out = C,
        idt_AB = idt_AB, idt_AC = idt_AC, idt_BC = idt_BC,
        mrca_node = u
      )
    }
  }
}

trios_df <- bind_rows(trios) %>%
  mutate(
    pass_range   = (idt_AB >= INGROUP_ID_MIN & idt_AB <= INGROUP_ID_MAX) &
      (idt_AC >= OUTGRP_ID_MIN  & idt_AC <= OUTGRP_ID_MAX)  &
      (idt_BC >= OUTGRP_ID_MIN  & idt_BC <= OUTGRP_ID_MAX)  &
      (pmin(idt_AC, idt_BC) < idt_AB),
    pass_balance = abs(idt_AC - idt_BC) <= BALANCE_DELTA,
    ratio        = parsimony_ratio(idt_AB, idt_AC, idt_BC)
  ) %>%
  filter(pass_range, pass_balance, !is.na(ratio)) %>%
  arrange(desc(ratio))

# --- 3) best trio per ingroup clade -----------------------------------------
best_per_clade <- trios_df %>%
  group_by(mrca_node) %>%
  slice_max(order_by = ratio, n = 1, with_ties = FALSE) %>%
  ungroup()
# Pick up the top ones if the values are equivalent.

# --- 4) decent clades among entire list -----------------------------------------
# knobs
RATIO_MIN <- 1.5   # keep all trios with ratio >= this criteria🐸
TOP_K     <- 3     # also cap to top K per ingroup pair (set Inf to disable)

trios_multi_per_pair <- trios_df %>%
  group_by(mrca_node, ing1, ing2) %>%
  arrange(desc(ratio), abs(idt_AC - idt_BC), desc(idt_AB)) %>%
  filter(ratio >= RATIO_MIN) %>%
  slice_head(n = TOP_K) %>%
  ungroup()

# --- 5) save outputs --------------------------------------------------------
readr::write_tsv(trios_df, "candidate_trios_all.tsv")
readr::write_tsv(best_per_clade, "candidate_trios_best_per_clade.tsv")

head(best_per_clade, 10)
