# Plan: Chi-square statistical tests for reviewer response

## Context

A reviewer requests statistical significance analyses that consider (a) deviations in empirical trinucleotide frequencies and (b) deviations of parsimony reconstruction from a Markov model of substitutions. Two chi-square tests are proposed:

1. **Spectrum comparison** — are two species' substitution spectra significantly different?
2. **Context dependence** — do flanking bases significantly affect the central-base substitution rate?

The tests should account for empirical trinucleotide frequencies (already captured in `totalRootNum` per TSV row) and address the Markov-model concern by noting the >80% identity criterion that ensures parsimony reliability.

---

## Statistical Design

### Key data structure (TSV files)
Each TSV (output of `trisbst_2TSVs.py`) has 96 rows, one per trinucleotide substitution type `X[Y>Z]W`:
- `mutNum`: observed substitution count
- `totalRootNum`: total sites of that trinucleotide type (= empirical frequency denominator)

The rate-normalized spectrum plotted in the thesis is `mutNum / totalRootNum`.

---

### Test 1 — `chi2_spectrum_cmp.py`: Compare two spectra

**Null hypothesis:** For every one of 96 trinucleotide types *i*, the substitution rate in species A equals that of species B.

**Method:** Chi-square test of homogeneity with **rate-adjusted expected counts**

For type *i*, pooled rate:
```
r_hat[i] = (mutNumA[i] + mutNumB[i]) / (rootNumA[i] + rootNumB[i])
```

Expected counts:
```
E[A, i] = rootNumA[i] × r_hat[i]
E[B, i] = rootNumB[i] × r_hat[i]
```

Chi-square statistic:
```
chi2 = Σ_i [(mutNumA[i] - E[A,i])² / E[A,i]  +  (mutNumB[i] - E[B,i])² / E[B,i]]
```

Degrees of freedom: 96 (= 2 species × 96 types − 96 pooled-rate estimates)

Why this is correct:
- Using `rootNum` as the exposure accounts for empirical trinucleotide frequencies (addresses reviewer point a)
- The pooled-rate null is equivalent to asking: "after correcting for trinucleotide composition, do the two spectra have the same shape?"
- The >80% identity criterion ensures parsimony reliability (addresses reviewer point b); `isParsimonious.py` provides per-trio parsimony-vs-non-parsimony ratios that can be cited

**Effect size:** Cramér's V = sqrt(chi2 / (N_total × 1)) where N_total = sum of all mutNum in both files. With millions of sites, p ≈ 0 even for small differences, so Cramér's V is the interpretively important number.

**CLI interface:**
```
python chi2_spectrum_cmp.py <tsv_file_A> <tsv_file_B> [--label-a NAME] [--label-b NAME]
```

**Output (stdout):**
```
chi2 = 12345.67   df = 96   p = 3.4e-123   V = 0.023
```
Optionally: per-type z-score to show which specific types drive the difference (useful for targeted biological interpretation).

---

### Test 2 — `chi2_context.py`: Context (flanking base) dependence

**Null hypothesis:** For a given central mutation type (e.g., C>A), the substitution rate is the same across all 16 flanking contexts.

**Method:** Chi-square goodness-of-fit for each of the 6 central mutation types

For central type *m*, pooled rate:
```
r_hat[m] = sum_i(mutNum[m,i]) / sum_i(rootNum[m,i])     (over 16 contexts i)
```

Expected counts:
```
E[m, i] = rootNum[m, i] × r_hat[m]
```

Chi-square statistic per central type:
```
chi2[m] = Σ_i (mutNum[m,i] - E[m,i])² / E[m,i]
```

Degrees of freedom: 15 (= 16 contexts − 1 estimated rate)

**Multiple testing correction:** 6 tests (one per central mutation type) → Bonferroni-corrected threshold p < 0.05/6 ≈ 0.0083. All 6 are expected to be highly significant; the interesting biological comparison is the chi2 magnitude (effect size) across types.

**Effect size:** η² = chi2 / (chi2 + N[m]) where N[m] = total sites for that central type.

**CLI interface:**
```
python chi2_context.py <tsv_file> [--label NAME]
```

**Output (stdout, one row per central mutation type):**
```
mutType   chi2       df   p           eta2
C>A       28431.2    15   <1e-300     0.041
C>G        4231.7    15   2.3e-045    0.009
...
```

---

## Files to Create

| File | Purpose |
|------|---------|
| `scripts/analysis/chi2_spectrum_cmp.py` | Test 1 script |
| `scripts/analysis/chi2_context.py` | Test 2 script |
| `scripts/analysis/test_chi2_spectrum_cmp.py` | Unit tests for Test 1 |
| `scripts/analysis/test_chi2_context.py` | Unit tests for Test 2 |
| `scripts/analysis/test/chi2_cmp_A.tsv` | Small fixture TSV for Test 1 |
| `scripts/analysis/test/chi2_cmp_B.tsv` | Small fixture TSV for Test 1 |
| `scripts/analysis/test/chi2_ctx.tsv` | Small fixture TSV for Test 2 |

No existing files need to be modified.

---

## Implementation Notes

- **Dependencies:** `scipy` (for `scipy.stats.chi2.sf`); `csv` from stdlib for TSV parsing
- **Style:** 4-space indent, snake_case, brief docstrings, argparse for CLI — consistent with existing analysis scripts
- **Edge cases to handle:** zero expected counts (skip type if E = 0); missing types in one file (treat as 0 mutNum)
- Both scripts share a common TSV-loading helper; extract to a small internal function (not a separate module — avoid over-engineering)
- Both scripts run independently with no dependency on MAF files or other pipeline steps

---

## Verification

```bash
cd scripts/analysis

# Run unit tests
python -m unittest test_chi2_spectrum_cmp test_chi2_context

# Smoke-test on real data (example paths)
python chi2_spectrum_cmp.py \
    ../data/ostMed_ostTau_ostLuc/ostMed2ostTau_one2one_20240423.tsv \
    ../data/ostMed_ostTau_ostLuc/ostMed2ostLuc_one2one_20240423.tsv

python chi2_context.py \
    ../data/ostMed_ostTau_ostLuc/ostMed2ostTau_one2one_20240423.tsv
```

Expected: chi-square values printed to stdout with finite p-values and Cramér's V / η² in the range 0.001–0.1 for real biological data.
