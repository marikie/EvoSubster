#!/usr/bin/env python3
"""Parse /usr/bin/time -v output files and produce a TSV summary table.

Usage:
    python parse_bench_times.py <time_dir> [--output bench_summary.tsv]

Expects files named: <DATE>_<LABEL>.time  (e.g. bench_20260412_agaBis_agaBit_agaSin.time)
Outputs one row per trio with wall-clock time, CPU%, memory, and exit status.
"""

import argparse
import os
import re
import sys


LINEAGE_MAP = {
    "agaBis", "armBor", "bolBar", "bolNan", "bolRex", "bolSem", "bolVar",
    "inoSue", "lacAka", "lacAme", "lacSan", "lecGlu", "lecObs", "lenEdo",
    "pleOst", "podMar", "podMin", "rusAbi", "strPac",
}

FUNGI_PREFIXES = {
    "agaBis", "armBor", "bolBar", "bolNan", "bolRex", "bolSem", "bolVar",
    "inoSue", "lacAka", "lacAme", "lacSan", "lecGlu", "lecObs", "lenEdo",
    "pleOst", "podMar", "podMin", "rusAbi", "strPac",
}

CNIDARIA_PREFIXES = {
    "acrAus", "acrDig", "acrMil", "casOrn", "cypSal", "dunAxi", "echHor",
    "hetMag", "monCap", "palCar", "palMiz", "pocGra", "porRus", "styPis",
}


def parse_time_file(path):
    """Extract metrics from a single /usr/bin/time -v output file."""
    fields = {
        "wall_clock": None,
        "user_sec": None,
        "system_sec": None,
        "cpu_percent": None,
        "max_rss_kb": None,
        "exit_status": None,
    }
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("Elapsed (wall clock) time"):
                    m = re.search(r":\s+(.+)$", line)
                    if m:
                        fields["wall_clock"] = m.group(1).strip()
                elif line.startswith("User time (seconds)"):
                    m = re.search(r":\s+([\d.]+)", line)
                    if m:
                        fields["user_sec"] = m.group(1)
                elif line.startswith("System time (seconds)"):
                    m = re.search(r":\s+([\d.]+)", line)
                    if m:
                        fields["system_sec"] = m.group(1)
                elif line.startswith("Percent of CPU this job got"):
                    m = re.search(r":\s+([\d]+)%", line)
                    if m:
                        fields["cpu_percent"] = m.group(1)
                elif line.startswith("Maximum resident set size"):
                    m = re.search(r":\s+([\d]+)", line)
                    if m:
                        fields["max_rss_kb"] = m.group(1)
                elif line.startswith("Exit status"):
                    m = re.search(r":\s+([\d]+)", line)
                    if m:
                        fields["exit_status"] = m.group(1)
    except OSError as e:
        print(f"Warning: could not read {path}: {e}", file=sys.stderr)
    return fields


def guess_lineage(label):
    prefix = label.split("_")[0]
    if prefix in FUNGI_PREFIXES:
        return "fungi"
    if prefix in CNIDARIA_PREFIXES:
        return "cnidaria"
    return "unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("time_dir", help="Directory containing .time files")
    parser.add_argument("--output", default=None,
                        help="Output TSV path (default: <time_dir>/bench_summary.tsv)")
    args = parser.parse_args()

    time_dir = args.time_dir
    if not os.path.isdir(time_dir):
        sys.exit(f"Error: {time_dir} is not a directory")

    out_path = args.output or os.path.join(time_dir, "bench_summary.tsv")

    # Collect and sort .time files
    time_files = sorted(
        f for f in os.listdir(time_dir) if f.endswith(".time")
    )
    if not time_files:
        sys.exit(f"No .time files found in {time_dir}")

    header = [
        "trio_label", "lineage",
        "wall_clock", "user_sec", "system_sec",
        "cpu_percent", "max_rss_kb", "exit_status",
    ]

    rows = []
    for fname in time_files:
        # Strip date prefix: <DATE>_<LABEL>.time  →  <LABEL>
        base = fname[:-5]  # remove .time
        parts = base.split("_", 1)
        label = parts[1] if len(parts) == 2 else base

        fields = parse_time_file(os.path.join(time_dir, fname))
        lineage = guess_lineage(label)
        rows.append([
            label, lineage,
            fields["wall_clock"] or "NA",
            fields["user_sec"] or "NA",
            fields["system_sec"] or "NA",
            fields["cpu_percent"] or "NA",
            fields["max_rss_kb"] or "NA",
            fields["exit_status"] or "NA",
        ])

    with open(out_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")

    print(f"Parsed {len(rows)} entries → {out_path}")


if __name__ == "__main__":
    main()
