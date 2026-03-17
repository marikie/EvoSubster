"""Compare two trinucleotide substitution spectra using a chi-square test of homogeneity.

Input:
    two TSV files (output of trisbst_2TSVs.py), one per species/dataset

Output:
    chi-square statistic, degrees of freedom, p-value, Cramér's V

Method: rate-adjusted expected counts.
For each trinucleotide type i:
    r_hat[i] = (mutNumA[i] + mutNumB[i]) / (rootNumA[i] + rootNumB[i])
    E[A,i]   = rootNumA[i] * r_hat[i]
    E[B,i]   = rootNumB[i] * r_hat[i]
    chi2    += (O_A - E_A)^2/E_A + (O_B - E_B)^2/E_B
df = number of trinucleotide types where both expected counts are nonzero

Using rootNum as exposure accounts for empirical trinucleotide frequencies.
The >80% pairwise identity criterion used in trio selection ensures that
parsimony-based substitution assignments are reliable.
"""

import argparse
import csv
import math

from scipy.stats import chi2 as chi2_dist


def load_tsv(path):
    """Return {mutType: (mutNum, totalRootNum)} from a trinucleotide TSV file."""
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            data[row["mutType"]] = (float(row["mutNum"]), float(row["totalRootNum"]))
    return data


def compute_chi2(data_a, data_b):
    """Compute chi-square statistic comparing two substitution spectra.

    Returns (chi2_stat, df, p_value, cramers_v, per_type).
    per_type: list of (mutType, chi2_contribution, sign)
        sign > 0 when rate in dataset A exceeds that in dataset B.
    """
    all_types = sorted(set(data_a) | set(data_b))
    chi2_stat = 0.0
    df = 0
    total_mut = 0.0
    per_type = []

    for mut_type in all_types:
        mut_a, root_a = data_a.get(mut_type, (0.0, 0.0))
        mut_b, root_b = data_b.get(mut_type, (0.0, 0.0))

        total_root = root_a + root_b
        if total_root == 0.0:
            continue

        r_hat = (mut_a + mut_b) / total_root
        exp_a = root_a * r_hat
        exp_b = root_b * r_hat

        # skip types where expected count is zero in either species
        if exp_a == 0.0 or exp_b == 0.0:
            continue

        contrib = (mut_a - exp_a) ** 2 / exp_a + (mut_b - exp_b) ** 2 / exp_b

        rate_a = mut_a / root_a if root_a > 0.0 else 0.0
        rate_b = mut_b / root_b if root_b > 0.0 else 0.0
        sign = 1 if rate_a >= rate_b else -1

        chi2_stat += contrib
        df += 1
        total_mut += mut_a + mut_b
        per_type.append((mut_type, contrib, sign))

    p_value = chi2_dist.sf(chi2_stat, df) if df > 0 else float("nan")
    cramers_v = math.sqrt(chi2_stat / total_mut) if total_mut > 0.0 else float("nan")

    return chi2_stat, df, p_value, cramers_v, per_type


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Chi-square test of homogeneity between two trinucleotide substitution "
            "spectra. Accounts for empirical trinucleotide frequencies via totalRootNum."
        )
    )
    parser.add_argument("tsv_a", help="TSV file for species/dataset A")
    parser.add_argument("tsv_b", help="TSV file for species/dataset B")
    parser.add_argument("--label-a", default="A", help="label for dataset A")
    parser.add_argument("--label-b", default="B", help="label for dataset B")
    parser.add_argument(
        "--per-type",
        action="store_true",
        help="print per-type chi2 contributions sorted by magnitude",
    )
    args = parser.parse_args()

    data_a = load_tsv(args.tsv_a)
    data_b = load_tsv(args.tsv_b)
    chi2_stat, df, p_value, cramers_v, per_type = compute_chi2(data_a, data_b)

    print(f"{args.label_a} vs {args.label_b}")
    print(f"chi2 = {chi2_stat:.4f}   df = {df}   p = {p_value:.4e}   V = {cramers_v:.4f}")

    if args.per_type:
        print()
        direction_label = f"direction ({args.label_a} higher = +)"
        print(f"{'mutType':<12}  {'chi2_contrib':>12}  {direction_label}")
        for mut_type, contrib, sign in sorted(per_type, key=lambda x: -x[1]):
            direction = "+" if sign > 0 else "-"
            print(f"{mut_type:<12}  {contrib:>12.4f}  {direction}")
