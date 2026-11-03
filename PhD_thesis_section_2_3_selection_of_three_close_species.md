# 2.3 Selection of three close species

To apply a parsimony-based framework, I selected triplets of closely related species in which pairwise sequence divergence is low. The outgroup was chosen to be close enough for parsimony to remain reliable, but sufficiently diverged from the other two species that the species-tree topology ((B, C), A) is unambiguous. This is important because mis-specifying the outgroup can substantially change the inferred number of substitutions assigned to each lineage.

To formalize this, consider a species tree with three species: A, B, and C, where A is the outgroup, and B and C form the ingroup.

Let $\mathrm{Id}_{AB}$, $\mathrm{Id}_{AC}$, and $\mathrm{Id}_{BC}$ denote the pairwise sequence identities between species A and B, A and C, and B and C, respectively. Here, sequence identity is defined as the proportion of aligned sites at which the two sequences share the same nucleotide, e.g.

$$
\mathrm{Id}_{AB} =
\frac{\#\{\text{aligned sites where A and B have the same base}\}}
{\#\{\text{aligned sites between A and B}\}}
$$

The corresponding mismatch rates are then

- between A and B: $x = 1 - \mathrm{Id}_{AB}$,
- between A and C: $y = 1 - \mathrm{Id}_{AC}$,
- between B and C: $z = 1 - \mathrm{Id}_{BC}$.

Let $p_A$, $p_B$, and $p_C$ denote the probabilities of mutation along the branches from the root to A, B, and C, respectively. Under a simple additive model of divergence, we can write

- $x = p_A + p_B$
- $y = p_A + p_C$
- $z = p_B + p_C$

which gives

- $p_A = \frac{x + y - z}{2}$
- $p_B = \frac{x + z - y}{2}$
- $p_C = \frac{y + z - x}{2}$.

For each ingroup species (B and C), I then compared the probability that a site pattern is explained by a single substitution (parsimony-consistent) versus two substitutions (non-parsimonious), under the assumption that any of the three non-ancestral nucleotides is equally likely (factor 1/3).

**Species B (ingroup)**

- Probability of a parsimony-consistent scenario (one substitution on branch B):

$$
p_{\mathrm{pars}}(B) = (1 - p_A)\frac{p_B}{3}(1 - p_C)
$$

- Probability of a non-parsimonious scenario (two substitutions on branches A and C):

$$
p_{\mathrm{nonpars}}(B) = \frac{p_A}{3}(1 - p_B)\frac{p_C}{3}
$$

Thus, the ratio of parsimony-consistent to non-parsimonious scenarios for species B is

$$
\mathrm{ratio}_B = 3 \cdot \frac{(1 - p_A)p_B(1 - p_C)}{p_A(1 - p_B)p_C}
$$

**Species C (ingroup)**

- Probability of a parsimony-consistent scenario (one substitution on branch C):

$$
p_{\mathrm{pars}}(C) = (1 - p_A)(1 - p_B)\frac{p_C}{3}
$$

- Probability of a non-parsimonious scenario (two substitutions on branches A and B):

$$
p_{\mathrm{nonpars}}(C) = \frac{p_A}{3}\frac{p_B}{3}(1 - p_C)
$$

The corresponding ratio for species C is

$$
\mathrm{ratio}_C = 3 \cdot \frac{(1 - p_A)(1 - p_B)p_C}{p_A p_B(1 - p_C)}
$$

These calculations were used as a rough check to confirm that, for the chosen species triplets, parsimony-consistent histories for the ingroup branches (B and C) are expected to be much more frequent than non-parsimonious two-hit histories.

In practice, I selected three species whose pairwise sequence identities were all above 80%. To ensure a reliable outgroup assignment, exclusion is driven by a two-step rule: if the genus pattern matches the target configuration (species A differs from both B and C while B and C share a genus), the trio is always retained, assuming that the original taxonomic classification reflects greater divergence. Otherwise, I examine substitution identities: when $\mathrm{Id}_{AB} < \mathrm{Id}_{BC}$ and $\mathrm{Id}_{AC} < \mathrm{Id}_{BC}$ hold, the trio is retained unless the genera form a two-vs-one configuration (species A shares a genus with only one ingroup species), in which case it is excluded.
