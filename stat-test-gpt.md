You: I'm under review for my journal paper. The pdf file loaded is my PhD thesis about the same research project as the paper I submitted. One reviewer said "there should be some support for statistical significance analyses. Note that the statistical modelling should consider the deviations in empirical frequencies and the deviations of parsimony reconstruction from a Markov model of substituions." What are "the deviations in empirical frequencies" and "the deviations of parsimony reconstruction from a Markov model of substitutions"? What kind of statistical test should I do?

Assistant: I think the reviewer is asking for **two different kinds of uncertainty** to be addressed, and they probably do **not** expect a single off-the-shelf test to cover both.

## 1) What they probably mean

### A. “Deviations in empirical frequencies”

This most likely means:

> your observed substitution frequencies are estimates from finite data, so they have sampling error / dispersion / uncertainty.

In your pipeline, a rate like

\[
\text{rate}\_{i}=\frac{\#\text{substitutions of type }i}{\#\text{ancestral contexts of type }i}
\]

is an **empirical frequency**.
The reviewer is saying: don’t just show the bar height — show whether that rate is meaningfully higher/lower than another rate **relative to its uncertainty**.

In your setting, this includes:

- finite counts,
- unequal numbers of ancestral contexts across substitution classes,
- overdispersion from genomic heterogeneity,
- local dependence because nearby sites are not fully independent,
- possible clustering in repeats / TE-rich regions.

So “deviations in empirical frequencies” = **how much the observed rates fluctuate around their underlying true rates**.

---

### B. “Deviations of parsimony reconstruction from a Markov model of substitutions”

This means:

> your substitution calls are inferred by parsimony, but the real process is more like a probabilistic continuous-time substitution process on a tree.

Parsimony says things like:

- if A=C≠B, infer a substitution on B.
- if A differs and B=C, exclude as ambiguous.
- if all differ, exclude.

But under a **Markov substitution model** (e.g. JC, HKY, GTR), the same observed leaf pattern can arise from:

- one substitution on branch B,
- or two substitutions on A and C,
- or other multi-hit histories.

So the reviewer is saying:

- parsimony is an approximation,
- quantify how much it may differ from a probabilistic model.

This is especially relevant for your trios near the lower identity end (~82–85%), where multiple hits are still uncommon but not negligible.

---

## 2) In plain language: what are they asking you to add?

Probably this:

1. **Statistical support for observed rate differences**
   e.g. “TTA>TCA is significantly enriched in Agaricomycetes”,
   “ACA>AAA is significantly enriched in cnidarians”,
   “CG>TG is elevated in sponges”.

2. **A sensitivity analysis showing that parsimony is not badly biasing the results**, or a correction based on a Markov model.

---

# 3) What test should you do?

## Short answer

For your paper, the most defensible combination would be:

### For empirical frequency differences:

Use a **count-based regression model with an offset for ancestral context counts**, ideally with overdispersion:

- **beta-binomial**
  or
- **negative binomial / quasi-Poisson GLM/GLMM with log(offset = context count)**

### For whole-spectrum differences:

Use:

- **PERMANOVA** on transformed spectra
  or
- a **global likelihood-ratio test** from a count model with context × lineage interaction

### For parsimony vs Markov:

Do a **sensitivity analysis** using either:

1. **maximum-likelihood ancestral reconstruction** under HKY/GTR, or
2. **parametric bootstrap simulations** under a fitted Markov model, then run your parsimony pipeline on the simulated data and measure bias.

---

# 4) What I would recommend specifically for your paper

You have **three kinds of claims** in the thesis/paper:

1. **specific context enrichments**
   - e.g. TTA>TCA and TTA>TGA in mushroom-forming fungi
   - ACA>AAA and ACG>AAG in cnidarians
   - CG>TG enrichment across lineages
   - unusual dinucleotide substitutions

2. **global spectrum differences among clades**
   - fungi vs cnidarians vs arthropods etc.
   - PCA separation

3. **validity of parsimony inference**

These should be handled separately.

---

## A. For specific substitution signatures: use a count model

For each substitution type \(i\) in genomic block \(w\) of species/trio \(s\):

- \(y\_{swi}\) = number of inferred substitutions of type \(i\)
- \(n\_{swi}\) = number of ancestral contexts/opportunities for type \(i\)

Then fit something like:

\[
y*{swi} \sim \text{NegBin}(\mu*{swi})
\]

\[
\log(\mu*{swi}) = \log(n*{swi}) + \alpha_i + \beta_i \cdot \text{group}\_s + u_s
\]

where:

- \(\log(n\_{swi})\) is an **offset**,
- \(\text{group}\_s\) is e.g. Agaricomycetes vs other fungi, or cnidarians vs non-cnidarians,
- \(u_s\) is a random effect for species/trio if you use a mixed model.

### Why this is good

It directly models:

- different numbers of opportunities,
- unequal counts,
- overdispersion,
- group-level enrichment.

### Practical interpretation

If \(\beta_i > 0\) and significant for TTA>TCA in Agaricomycetes, you can say:

> the rate of TTA>TCA is significantly elevated in Agaricomycetes after accounting for the number of TTA ancestral contexts.

---

## B. Use genomic blocks, not one genome-wide number only

This is important.

If you use one genome-wide count per species, the p-values can become misleading because:

- you have no within-genome variance estimate,
- sites are not independent,
- overlap among k-mers adds dependence.

### Better:

Split the alignment into **non-overlapping blocks** (e.g. 100 kb, 500 kb, or 1 Mb), and compute:

- substitution counts per block,
- ancestral context counts per block.

Then use those blocks as replicate observations.

This gives you:

- confidence intervals,
- better error estimation,
- robustness to local clustering.

If you want a simpler nonparametric approach, use a **block bootstrap**:

- resample blocks with replacement,
- recompute rates / differences,
- get 95% CIs.

That is much better than naive binomial SEs.

---

## C. For overall spectrum differences: global test

If you want to support PCA/clustering with a p-value, use:

### Option 1: PERMANOVA

On:

- log-transformed or centered-log-ratio transformed spectra
- distance = Euclidean on CLR, or Bray–Curtis / Aitchison

Test:

- do spectra differ among major groups (e.g. phyla, classes)?

Then also test **homogeneity of dispersion** (e.g. `betadisper` in vegan), because PERMANOVA can be confounded by dispersion differences.

### Option 2: Global likelihood-ratio test from the count model

Fit:

- model without group × context interaction
- model with group × context interaction

Then compare by LRT.

Interpretation:

- if interaction is significant, the overall substitution spectrum differs among groups.

This is arguably more aligned with your count data than PERMANOVA.

---

# 5) What about the “parsimony vs Markov” part?

This is the part many reviewers care about conceptually.

## Best practical solution: Markov-model sensitivity analysis

### Step 1: Fit a simple substitution model to each trio

Use neutral/noncoding aligned sites and fit, for each trio:

- **HKY85** or **GTR**
- with branch lengths on the fixed topology \(((B,C),A)\)

You do **not** need a context-dependent 96-state model just to satisfy this reviewer.
A standard nucleotide CTMC is probably enough as a sensitivity analysis.

### Step 2: Quantify parsimony bias in one of two ways

---

## Option A: ML ancestral reconstruction

For each site pattern, under the fitted Markov model, compute the posterior probability that the substitution occurred on branch B or C.

Then instead of hard parsimony counts, you get **expected counts**:

- parsimony count = 1 or 0
- Markov count = posterior-weighted expected contribution

Compare:

- parsimony-inferred counts vs Markov-model expected counts
- per trio and per substitution type

If results are very similar, that strongly helps your paper.

---

## Option B: Parametric bootstrap simulation

This is often the easiest to explain to reviewers.

For each trio:

1. fit HKY/GTR + branch lengths,
2. simulate alignments of the same length under that model,
3. run your existing parsimony pipeline on those simulated alignments,
4. compare:
   - true simulated branch-specific substitutions
   - parsimony-inferred substitutions

This gives you:

- bias,
- variance,
- error rates of parsimony under realistic branch lengths.

Then you can say something like:

> Across simulated alignments under trio-specific HKY models, parsimony reconstruction closely approximated the true branch-specific substitution counts over the divergence range analyzed, with low context-specific bias for the principal signatures discussed in the manuscript.

That directly answers the reviewer.

---

# 6) How to connect this to your thesis specifically

Your thesis already has a **very useful starting point** in section 2.3:

- you derive trio-specific branch probabilities from pairwise identities,
- and compute ratios of parsimony-consistent vs non-parsimonious histories.

That is already the beginning of a **parsimony error analysis**.

So you can expand that into something more explicitly statistical.

## Minimal analytical version

For each trio, compute the expected probability that a site pattern assigned by parsimony to branch B was actually generated by a non-parsimonious history under a simple Markov-style model.

Using your notation, if you already have:

\[
\text{ratio}\_B = \frac{P(\text{parsimonious one-hit on B})}{P(\text{non-parsimonious two-hit alternative})}
\]

then an approximate misassignment probability is:

\[
P(\text{parsimony wrong} \mid \text{pattern assigned to B})
\approx \frac{1}{1+\text{ratio}\_B}
\]

and similarly for C.

This is not as strong as full HKY/GTR simulation, but it is much better than nothing, and it is closely aligned with your current framework.

---

# 7) What I think is the best “revision-friendly” plan

If you want something strong but still realistic for a revision, I would do this:

## Recommended package of analyses

### 1. For key highlighted signatures

Use **block-based negative binomial or beta-binomial models** with ancestral context count as offset/denominator.

Examples:

- TTA>TCA in Agaricomycetes vs other fungi
- TTA>TGA in Agaricomycetes vs other fungi
- ACA>AAA and ACG>AAG in cnidarians vs non-cnidarians
- CG>TG in sponges vs other metazoans
- TT>AA / AT>TA / TA>AT in Glomeromycotina vs other fungi

Report:

- rate ratio,
- 95% CI,
- p-value,
- FDR-adjusted p-value.

### 2. For global spectrum differences

Use either:

- **PERMANOVA** on transformed spectra, or
- a **global context × lineage interaction test** in a count model.

### 3. For parsimony validity

Do **trio-specific HKY/GTR simulation-based sensitivity analysis**.

Report:

- correlation between parsimony and model-based counts,
- average absolute deviation or relative bias,
- maybe a supplementary figure/table.

This would satisfy the reviewer much better than just adding a few chi-square tests.

---

# 8) What not to do

I would **not** recommend:

- simple t-tests on normalized bar heights,
- plain chi-square tests on raw substitution counts,
- Fisher tests ignoring unequal context opportunities,
- p-values based only on total genome-wide counts without block structure,
- saying “all pairwise identities are above 80%, therefore parsimony is fine” without quantitative support.

Those are weaker than what the reviewer is explicitly asking for.

---

# 9) A concrete model you can write in Methods

Here is a clean formulation.

## For individual substitution types

For substitution type \(i\), block \(w\), species \(s\):

- \(y\_{swi}\): inferred substitutions
- \(n\_{swi}\): ancestral context counts

Model:

\[
y*{swi} \sim \text{NegBin}(\mu*{swi}, \theta_i)
\]

\[
\log(\mu*{swi}) = \log(n*{swi}) + \alpha_i + \beta_i G_s + u_s
\]

where:

- \(G_s\) is lineage/group,
- \(u_s\) is a species/trio random effect.

Test:
\[
H_0:\beta_i = 0
\]

Correct across many substitution types with **Benjamini–Hochberg FDR**.

---

## For overall spectra

Fit a model with context × group interaction:

\[
\log(\mu*{swi}) = \log(n*{swi}) + \alpha*i + \gamma*{G*s} + \delta*{i,G_s} + u_s
\]

Then test:

\[
H*0: \delta*{i,G_s}=0 \text{ for all } i
\]

using a likelihood-ratio test.

That tests whether the **shape of the spectrum** differs among groups.

---

# 10) If you want the simplest acceptable addition

If time is tight, do this minimum:

### Minimum acceptable statistical response

1. **Block bootstrap** 95% CIs for major substitution rates and key lineage contrasts.
2. **Negative binomial / beta-binomial tests** for the specific highlighted signatures.
3. **A parsimony sensitivity analysis under a fitted HKY/GTR model**, ideally by simulation.

That is probably enough for the reviewer.

---

# 11) Suggested wording for your response to the reviewer

You could say something like:

> We agree that statistical support should reflect both sampling variability in empirical substitution frequencies and the approximation introduced by parsimony-based inference. We therefore added two analyses. First, we estimated uncertainty in context-specific substitution rates using block-wise count data and generalized count models with the number of ancestral contexts as exposure, reporting confidence intervals and FDR-adjusted significance for key lineage-specific signatures. Second, because parsimony is an approximation to a continuous-time Markov substitution process, we performed trio-specific sensitivity analyses under fitted HKY/GTR models, using probabilistic reconstruction and/or parametric bootstrap simulations to quantify the deviation of parsimony-based counts from model-based expectations. The main qualitative conclusions were robust.

---

# 12) My bottom-line recommendation

For **your project**, the best answer is:

- **Do not rely on one single test.**
- Use **count-based models with context-count offsets** for significance of specific signatures.
- Use **PERMANOVA or a global interaction test** for whole-spectrum differences.
- Use **HKY/GTR-based simulation or ML reconstruction** to show that parsimony bias is small enough over your trio divergence range.

If you want, I can help you design:

1. the exact statistical analysis section for the paper,
2. an R workflow (`glmmTMB` / `vegan` / bootstrap), or
3. a reviewer response paragraph tailored to your manuscript.
