#!/usr/bin/env python3

"""
Apply the thesis (section 2.3) trio selection rule.

A trio is one outgroup plus two ingroup species, and the pipeline assumes the
topology ((org2, org3), org1) -- org1 is the outgroup, org2 and org3 are the
ingroup.  A trio is retained only when both the pairwise sequence identities and
the genus names are consistent with that topology:

  * every pairwise identity must exceed the threshold (default 80%);
  * if the genera form the target configuration (the outgroup differs from both
    ingroup species, which share a genus) the trio is always retained;
  * otherwise the identities must order as idt_12 < idt_23 and idt_13 < idt_23
    (the ingroup pair is the most similar), and the trio is still excluded when
    the genera split two-vs-one, i.e. the outgroup shares a genus with exactly
    one ingroup species -- taxonomy is then asserting a topology that
    contradicts the outgroup assignment.

This module is the single source of that rule.  src/report/collect_run_summary.py
imports it for per-run reporting; trio_selection.R shells out to the batch CLI
below to screen candidate trios enumerated from a phylogenetic tree.

Batch CLI:

    python3 src/select/trio_filter.py --trios trios.tsv --out verdict.tsv

Input TSV requires the columns idt_12, idt_13, idt_23 (percent identities) and
genus_1, genus_2, genus_3 (organism or genus names, slot 1 = outgroup).  Any
other columns are passed through untouched.  The output repeats the input
columns and appends every key of filter_details plus filter_issues.  Booleans are
written as TRUE/FALSE and unknowns as NA so that R reads them natively.
"""

import argparse
import csv
import re
import sys
from typing import Dict, List, Optional, Tuple

# Slot 1 is the outgroup; slots 2 and 3 are the ingroup.  The thesis uses the
# opposite lettering (A = outgroup, B/C = ingroup), so its formulas cannot be
# transcribed literally.
SLOTS = ("org1", "org2", "org3")
IDENTITY_KEYS = ("idt_12", "idt_13", "idt_23")
GENUS_COLUMNS = ("genus_1", "genus_2", "genus_3")

FILTER_DETAIL_KEYS = (
    "genus_condition_met",
    "genus_pattern_same_all",
    "genus_pattern_all_different",
    "genus_pattern_two_vs_one",
    "identity_condition_met",
    "idt_threshold_condition",
    "idt_threshold_value",
    "filter_flag",
    "excluded",
)


def extract_genus_label(metadata_entry: dict) -> Optional[str]:
    candidates = [
        metadata_entry.get("raw_organism_name"),
        metadata_entry.get("ncbi_full_name"),
        metadata_entry.get("directory_name"),
        metadata_entry.get("short_name"),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        sanitized = candidate.replace("_", " ").strip()
        if not sanitized:
            continue
        token = sanitized.split()[0]
        if token:
            return token
    return None


def evaluate_genus_condition(slot_map: Dict[str, dict]) -> Tuple[Optional[bool], List[str]]:
    patterns, issues = evaluate_genus_patterns(slot_map)
    if patterns is None:
        return None, issues
    return patterns["pattern1"], issues


def evaluate_genus_patterns(
    slot_map: Dict[str, dict]
) -> Tuple[Optional[Dict[str, bool]], List[str]]:
    issues: List[str] = []
    required_slots = ("org1", "org2", "org3")
    genus_values: Dict[str, str] = {}

    for slot in required_slots:
        entry = slot_map.get(slot)
        if not entry:
            issues.append(f"{slot}: Missing metadata entry for genus evaluation.")
            return None, issues
        genus_label = extract_genus_label(entry)
        if not genus_label:
            issues.append(f"{slot}: Unable to derive genus from metadata.")
            return None, issues
        genus_values[slot] = genus_label.lower()

    pattern1 = (
        genus_values["org1"] != genus_values["org2"]
        and genus_values["org1"] != genus_values["org3"]
        and genus_values["org2"] == genus_values["org3"]
    )
    pattern3 = (
        genus_values["org1"] == genus_values["org2"]
        and genus_values["org2"] == genus_values["org3"]
    )
    pattern4 = (
        genus_values["org1"] != genus_values["org2"]
        and genus_values["org1"] != genus_values["org3"]
        and genus_values["org2"] != genus_values["org3"]
    )
    pattern5 = (
        (
            genus_values["org1"] != genus_values["org2"]
            and genus_values["org1"] == genus_values["org3"]
            and genus_values["org2"] != genus_values["org3"]
        )
        or (
            genus_values["org1"] == genus_values["org2"]
            and genus_values["org1"] != genus_values["org3"]
            and genus_values["org2"] != genus_values["org3"]
        )
    )

    patterns = {
        "pattern1": pattern1,
        "pattern3": pattern3,
        "pattern4": pattern4,
        "pattern5": pattern5,
    }
    return patterns, issues


def parse_identity_percentage(value: str) -> Optional[float]:
    if not value:
        return None
    sanitized = value.replace("%", " ")
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", sanitized)
    if not match:
        return None
    try:
        return float(match.group())
    except ValueError:
        return None


def evaluate_identity_condition(
    identity_values: Dict[str, str]
) -> Tuple[Optional[bool], List[str], Optional[Dict[str, float]]]:
    issues: List[str] = []
    numeric_values: Dict[str, float] = {}
    required_keys = ("idt_12", "idt_13", "idt_23")

    for key in required_keys:
        raw_value = identity_values.get(key)
        if raw_value is None:
            issues.append(f"{key}: Identity value missing; cannot evaluate condition.")
            return None, issues, None
        parsed = parse_identity_percentage(raw_value)
        if parsed is None:
            issues.append(f"{key}: Unable to parse numeric identity from '{raw_value}'.")
            return None, issues, None
        numeric_values[key] = parsed

    condition_met = numeric_values["idt_12"] < numeric_values["idt_23"] and numeric_values["idt_13"] < numeric_values["idt_23"]
    return condition_met, issues, numeric_values


def evaluate_filter_status(
    slot_map: Dict[str, dict],
    identity_values: Dict[str, str],
    idt_threshold: float,
) -> Tuple[Dict[str, Optional[bool]], List[str]]:
    issues: List[str] = []
    genus_patterns, genus_issues = evaluate_genus_patterns(slot_map)
    identity_condition, identity_issues, numeric_identity_values = evaluate_identity_condition(
        identity_values
    )
    issues.extend(genus_issues)
    issues.extend(identity_issues)

    genus_condition = genus_patterns["pattern1"] if genus_patterns else None
    genus_pattern_same_all = genus_patterns["pattern3"] if genus_patterns else None
    genus_pattern_all_different = genus_patterns["pattern4"] if genus_patterns else None
    genus_pattern_two_vs_one = genus_patterns["pattern5"] if genus_patterns else None

    threshold_condition: Optional[bool] = None
    if numeric_identity_values is not None:
        threshold_condition = all(value > idt_threshold for value in numeric_identity_values.values())

    filter_flag: Optional[bool] = None

    # Threshold gate: if below threshold, filter out before other rules
    if threshold_condition is False:
        filter_flag = True
    else:
        if genus_condition is True:
            filter_flag = False
        elif identity_condition is True:
            if genus_pattern_same_all is True:
                filter_flag = False
            elif genus_pattern_all_different is True:
                filter_flag = False
            elif genus_pattern_two_vs_one is True:
                filter_flag = True
        elif identity_condition is False:
            filter_flag = True

    excluded = filter_flag is True
    filter_details = {
        "genus_condition_met": genus_condition,
        "genus_pattern_same_all": genus_pattern_same_all,
        "genus_pattern_all_different": genus_pattern_all_different,
        "genus_pattern_two_vs_one": genus_pattern_two_vs_one,
        "identity_condition_met": identity_condition,
        "idt_threshold_condition": threshold_condition,
        "idt_threshold_value": idt_threshold,
        "filter_flag": filter_flag,
        "excluded": excluded,
    }
    return filter_details, issues


# --- batch CLI ---------------------------------------------------------------


def format_cell(value) -> str:
    """Render a filter_details value so that R reads it natively."""
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return str(value)


def screen_row(row: Dict[str, str], idt_threshold: float) -> Dict[str, str]:
    slot_map = {
        slot: {"raw_organism_name": row.get(column, "")}
        for slot, column in zip(SLOTS, GENUS_COLUMNS)
    }
    identity_values = {key: row.get(key) for key in IDENTITY_KEYS}

    filter_details, issues = evaluate_filter_status(slot_map, identity_values, idt_threshold)

    verdict = {key: format_cell(filter_details[key]) for key in FILTER_DETAIL_KEYS}
    verdict["filter_issues"] = "; ".join(issues)
    return verdict


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply the thesis section 2.3 trio selection rule to a table of candidate trios."
    )
    parser.add_argument(
        "--trios",
        required=True,
        help=(
            "Input TSV. Requires columns idt_12, idt_13, idt_23, genus_1, genus_2, genus_3 "
            "(slot 1 = outgroup). Other columns are passed through."
        ),
    )
    parser.add_argument(
        "--idt-threshold",
        type=float,
        default=80.0,
        help="Every pairwise identity must exceed this to be retained (default: 80).",
    )
    parser.add_argument(
        "--out",
        help="Output TSV path. Defaults to stdout.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    with open(args.trios, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            print(f"error: {args.trios} has no header row.", file=sys.stderr)
            return 1

        missing = [
            column
            for column in IDENTITY_KEYS + GENUS_COLUMNS
            if column not in reader.fieldnames
        ]
        if missing:
            print(
                f"error: {args.trios} is missing required column(s): {', '.join(missing)}",
                file=sys.stderr,
            )
            return 1

        rows = list(reader)

    out_columns = list(reader.fieldnames) + list(FILTER_DETAIL_KEYS) + ["filter_issues"]
    handle = open(args.out, "w", newline="", encoding="utf-8") if args.out else sys.stdout
    try:
        writer = csv.DictWriter(
            handle, fieldnames=out_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            row.update(screen_row(row, args.idt_threshold))
            writer.writerow(row)
    finally:
        if args.out:
            handle.close()

    retained = sum(1 for row in rows if row["excluded"] == "FALSE")
    print(
        f"{len(rows)} trios screened at idt-threshold {args.idt_threshold}: "
        f"{retained} retained, {len(rows) - retained} excluded.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
