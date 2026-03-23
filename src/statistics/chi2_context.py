"""Test whether flanking bases significantly affect the central-base substitution rate.

Input:
    one TSV file (output of trisbst_2TSVs.py)

Output:
    for each of the 6 central mutation types (C>A, C>G, C>T, T>A, T>C, T>G):
        chi-square statistic, degrees of freedom, p-value, eta-squared

Method: chi-square goodness-of-fit against a uniform-rate null hypothesis.
For each central mutation type m (e.g. C>A) and its flanking contexts i:
    r_hat[m] = sum_i(mutNum[m,i]) / sum_i(rootNum[m,i])
    E[m,i]   = rootNum[m,i] * r_hat[m]
    chi2[m]  = sum_i (mutNum[m,i] - E[m,i])^2 / E[m,i]
    df       = (number of valid contexts) - 1

Using rootNum as exposure accounts for empirical trinucleotide frequencies.
Bonferroni correction applies across the 6 central mutation types (threshold p < 0.05/6).
"""

import argparse
import csv
import math

from scipy.stats import chi2 as chi2_dist

BONFERRONI_N = 6  # number of central mutation types


def load_tsv(path):
    """Return {mutType: (mutNum, totalRootNum)} from a trinucleotide TSV file."""
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            data[row["mutType"]] = (float(row["mutNum"]), float(row["totalRootNum"]))
    return data


def central_mut(mut_type):
    """Extract central mutation type from 'X[Y>Z]W' format. Returns 'Y>Z'."""
    return mut_type[2:5]


def group_by_central(data):
    """Group TSV data by central mutation type.

    Returns {central_mut: {mutType: (mutNum, totalRootNum)}}.
    """
    groups = {}
    for mut_type, counts in data.items():
        cm = central_mut(mut_type)
        if cm not in groups:
            groups[cm] = {}
        groups[cm][mut_type] = counts
    return groups


def compute_chi2_context(contexts):
    """Compute chi-square goodness-of-fit for one central mutation type.

    contexts: {mutType: (mutNum, totalRootNum)}
    Returns (chi2_stat, df, p_value, eta2).
    """
    total_mut = sum(m for m, _ in contexts.values())
    total_root = sum(r for _, r in contexts.values())

    if total_root == 0.0 or total_mut == 0.0:
        return float("nan"), 0, float("nan"), float("nan")

    r_hat = total_mut / total_root

    chi2_stat = 0.0
    valid = 0
    for mut_num, root_num in contexts.values():
        if root_num == 0.0:
            continue
        exp = root_num * r_hat
        if exp == 0.0:
            continue
        chi2_stat += (mut_num - exp) ** 2 / exp
        valid += 1

    df = max(valid - 1, 0)
    p_value = chi2_dist.sf(chi2_stat, df) if df > 0 else float("nan")
    eta2 = chi2_stat / (chi2_stat + total_mut) if (chi2_stat + total_mut) > 0.0 else float("nan")

    return chi2_stat, df, p_value, eta2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Chi-square goodness-of-fit test for context dependence of substitution "
            "rates. Tests whether flanking bases significantly affect the substitution "
            "rate for each central mutation type."
        )
    )
    parser.add_argument("tsv", help="TSV file (output of trisbst_2TSVs.py)")
    parser.add_argument("--label", default="", help="dataset label for display")
    args = parser.parse_args()

    data = load_tsv(args.tsv)
    groups = group_by_central(data)

    if args.label:
        print(args.label)

    threshold = 0.05 / BONFERRONI_N
    header = (
        f"{'mutType':<8}  {'chi2':>12}  {'df':>4}  {'p':>12}  {'eta2':>8}  sig*"
    )
    print(header)
    print("-" * len(header))

    for cm in sorted(groups):
        chi2_stat, df, p_value, eta2 = compute_chi2_context(groups[cm])
        sig = "*" if (not math.isnan(p_value) and p_value < threshold) else ""
        p_str = f"{p_value:.4e}" if not math.isnan(p_value) else "nan"
        eta2_str = f"{eta2:.4f}" if not math.isnan(eta2) else "nan"
        print(
            f"{cm:<8}  {chi2_stat:>12.4f}  {df:>4}  {p_str:>12}  {eta2_str:>8}  {sig}"
        )

    print(f"\n* Bonferroni threshold: p < {threshold:.4f} (n = {BONFERRONI_N} tests)")
