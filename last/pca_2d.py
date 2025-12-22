#!/usr/bin/env python3

"""
Traverse dataset directories produced by trisbst workflows, gather substitution TSV files
from each latest run, and generate 2D PCA clustering plots along with CSV summaries.

Usage:
    python scripts/last/pca_2d.py /path/to/data/cnidaria --tsv-pattern "*_maflinked.tsv"
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text  # type: ignore
from matplotlib import rcParams
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

rcParams["font.family"] = "DejaVu Sans"

EPS = 1e-12
DEFAULT_TSV_PATTERN = "*_maflinked_ncds.tsv"
IDENTITY_PREFIX = "# substitution percent identity:"
CLASSIFICATION_ORDER = [
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Walk dataset directories under the input root, load substitution TSV files from each latest run, "
            "and run 2D PCA clustering on the resulting logRatio vectors."
    )
    )
    parser.add_argument(
        "input_dir",
        help="Root directory containing dataset subdirectories (e.g., data/cnidaria).",
    )
    parser.add_argument(
        "--tsv-pattern",
        default=DEFAULT_TSV_PATTERN,
        help=f"Glob pattern to select TSV files under each run directory (default: {DEFAULT_TSV_PATTERN}).",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory where PCA outputs will be written. Defaults to <input_dir>/pca.",
    )
    parser.add_argument(
        "--max-clusters",
        type=int,
        default=10,
        help="Maximum number of clusters to evaluate when selecting the optimal k.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed used for clustering.",
    )
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        parser.print_help(sys.stderr)
        print(
            "\nExample:\n  python scripts/last/pca_2d.py "
            "/home/mrk/sbst/data/cnidaria --tsv-pattern '*_maflinked.tsv' --output-dir /tmp/cnidaria_pca",
            file=sys.stderr,
        )
        raise SystemExit(1)
    return parser.parse_args(argv)


def parse_filename_parts(tsv_path: Path) -> Tuple[str, str, str]:
    name_parts = tsv_path.stem.split("_")
    if len(name_parts) < 3:
        raise ValueError(f"Unexpected filename format for TSV file: {tsv_path.name}")
    accession = "_".join(name_parts[:2])
    short_name = name_parts[2]
    date = name_parts[3] if len(name_parts) > 3 else ""
    return accession, short_name, date


def warn(message: str) -> None:
    print(f"Warning: {message}", file=sys.stderr)


def find_latest_run_dir(dataset_dir: Path) -> Optional[Path]:
    run_dirs = [
        entry for entry in dataset_dir.iterdir() if entry.is_dir() and entry.name.isdigit()
    ]
    if not run_dirs:
        return None
    return max(run_dirs, key=lambda path: path.name)


def resolve_path(value: str, base_dir: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def load_manifest_lookup(
    manifest_path: Path,
) -> Tuple[dict, Dict[Tuple[str, str], Dict[str, Path]], Path, Dict[str, dict]]:
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    metadata_dir_value = manifest.get("metadata_dir")
    if not metadata_dir_value:
        raise ValueError(f"'metadata_dir' is missing in manifest {manifest_path}")

    metadata_dir = resolve_path(metadata_dir_value, manifest_path.parent)
    if not metadata_dir.exists():
        raise FileNotFoundError(f"metadata_dir does not exist: {metadata_dir}")

    lookup: Dict[Tuple[str, str], Dict[str, Path]] = {}
    for organism in manifest.get("organisms", []):
        role = (organism.get("role") or "").lower()
        if role == "outgroup":
            continue
        short = organism.get("short_name")
        accession = organism.get("accession")
        metadata_json = organism.get("metadata_json")

        if not short or not accession or not metadata_json:
            raise ValueError(f"Incomplete organism entry in manifest {manifest_path}: {organism}")

        metadata_path = resolve_path(metadata_json, manifest_path.parent)
        taxonomy_path = metadata_dir / f"taxonomy_{short}_{accession}.json"
        lookup[(accession, short)] = {
            "metadata_path": metadata_path,
            "taxonomy_path": taxonomy_path,
        }

    if not lookup:
        raise ValueError(f"No ingroup organisms found in manifest {manifest_path}")

    slot_map = build_manifest_slot_map(manifest)
    return manifest, lookup, metadata_dir, slot_map


def load_json_cached(path: Path, cache: Dict[Path, dict]) -> dict:
    if path not in cache:
        if not path.exists():
            raise FileNotFoundError(f"Required JSON file not found: {path}")
        with path.open("r", encoding="utf-8") as handle:
            cache[path] = json.load(handle)
    return cache[path]


def extract_tax_id(metadata: dict) -> Optional[str]:
    reports = metadata.get("reports")
    if isinstance(reports, list) and reports:
        first = reports[0]
        organism = first.get("organism", {})
        tax_id = organism.get("tax_id")
        if tax_id:
            return str(tax_id)
        biosample = first.get("biosample", {}).get("description", {}).get("organism", {})
        tax_id = biosample.get("tax_id")
        if tax_id:
            return str(tax_id)
    return None


def ensure_taxonomy_json(taxonomy_path: Path, metadata: dict, short_name: str, accession: str) -> bool:
    if taxonomy_path.exists() and taxonomy_path.stat().st_size > 0:
        return True

    tax_id = extract_tax_id(metadata)
    if not tax_id:
        warn(f"{short_name}: tax_id not found in metadata; cannot fetch taxonomy JSON.")
        return False

    taxonomy_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        result = subprocess.run(
            ["datasets", "summary", "taxonomy", "taxon", tax_id],
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        warn("datasets CLI not found; skipping taxonomy download.")
        return False
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else "unknown error"
        warn(f"{short_name}: Failed to download taxonomy for tax_id {tax_id}: {stderr}")
        return False

    tmp_file: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False) as handle:
            handle.write(result.stdout)
            tmp_file = Path(handle.name)
        shutil.move(str(tmp_file), taxonomy_path)
        return True
    except OSError as exc:
        warn(f"{short_name}: Failed to write taxonomy JSON for tax_id {tax_id}: {exc}")
        if tmp_file and tmp_file.exists():
            tmp_file.unlink(missing_ok=True)
        return False


def extract_classification_map(taxonomy: dict, identifier: str) -> Dict[str, str]:
    reports = taxonomy.get("reports", [])
    if not reports:
        warn(f"taxonomy JSON missing 'reports' for {identifier}")
        return {}

    classification = reports[0].get("taxonomy", {}).get("classification", {})
    labels: Dict[str, str] = {}
    for level, entry in classification.items():
        if not isinstance(entry, dict):
            continue
        name = entry.get("name")
        if name:
            labels[level.lower()] = name

    if not labels:
        warn(f"taxonomy classification map missing for {identifier}")
    return labels


def extract_genus(metadata: dict, short_name: str, accession: str) -> Tuple[str, str]:
    reports = metadata.get("reports", [])
    if reports:
        organism = reports[0].get("organism", {})
        organism_name = organism.get("organism_name")
        if organism_name:
            genus = organism_name.split()[0]
            return genus, organism_name

    fallback_label = accession
    return "Unknown", fallback_label


def build_manifest_slot_map(manifest: dict) -> Dict[str, dict]:
    slot_map: Dict[str, dict] = {}
    for entry in manifest.get("organisms", []):
        slot = entry.get("slot")
        if slot:
            slot_map[slot] = entry
    return slot_map


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

    condition_met = (
        genus_values["org1"] != genus_values["org2"]
        and genus_values["org1"] != genus_values["org3"]
        and genus_values["org2"] == genus_values["org3"]
    )
    return condition_met, issues


def select_first(run_dir: Path, pattern: str) -> Optional[Path]:
    matches = sorted(run_dir.glob(pattern))
    if not matches:
        return None
    return matches[0]


def extract_identity_line(train_path: Path) -> Optional[str]:
    lines = [
        line.strip()
        for line in train_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(IDENTITY_PREFIX)
    ]
    if len(lines) >= 2:
        return lines[-2]
    if lines:
        return lines[0]
    return None


def extract_identity_value(identity_line: str) -> str:
    if ":" in identity_line:
        _, remainder = identity_line.split(":", 1)
        return remainder.strip()
    return identity_line.strip()


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


def collect_identity_strings(run_dir: Path, slot_map: Dict[str, dict]) -> Tuple[Dict[str, str], List[str]]:
    identity: Dict[str, str] = {}
    issues: List[str] = []

    def short_name_for(slot: str) -> Optional[str]:
        entry = slot_map.get(slot)
        return entry.get("short_name") if entry else None

    pair_specs = [
        ("idt_12", "org1", "org2"),
        ("idt_13", "org1", "org3"),
        ("idt_23", "org2", "org3"),
    ]

    for key, slot_a, slot_b in pair_specs:
        short_a = short_name_for(slot_a)
        short_b = short_name_for(slot_b)
        if not short_a or not short_b:
            issues.append(f"{key}: Missing short_name for {slot_a} or {slot_b}.")
            continue

        pattern = f"*{short_a}2{short_b}_*.train"
        train_path = select_first(run_dir, pattern)
        if not train_path:
            issues.append(f"{key}: No train file matching '{pattern}' found.")
            continue

        identity_line = extract_identity_line(train_path)
        if not identity_line:
            issues.append(f"{key}: '{IDENTITY_PREFIX}' line missing in {train_path.name}.")
            continue

        raw_value = extract_identity_value(identity_line)
        value_with_unit = raw_value if raw_value.endswith("%") else f"{raw_value} %"
        identity[key] = value_with_unit

    return identity, issues


def evaluate_identity_condition(identity_values: Dict[str, str]) -> Tuple[Optional[bool], List[str]]:
    if not identity_values:
        return None, []

    issues: List[str] = []
    numeric_values: Dict[str, float] = {}
    required_keys = ("idt_12", "idt_13", "idt_23")

    for key in required_keys:
        raw_value = identity_values.get(key)
        if raw_value is None:
            issues.append(f"{key}: Identity value missing; cannot evaluate condition.")
            return None, issues
        parsed = parse_identity_percentage(raw_value)
        if parsed is None:
            issues.append(f"{key}: Unable to parse numeric identity from '{raw_value}'.")
            return None, issues
        numeric_values[key] = parsed

    condition_met = numeric_values["idt_12"] < numeric_values["idt_23"] and numeric_values["idt_13"] < numeric_values["idt_23"]
    return condition_met, issues


def evaluate_dataset_filter(slot_map: Dict[str, dict], run_dir: Optional[Path]) -> Tuple[Dict[str, Optional[bool]], List[str]]:
    issues: List[str] = []
    genus_condition, genus_issues = evaluate_genus_condition(slot_map)
    issues.extend(genus_issues)

    identity_values: Dict[str, str] = {}
    if run_dir and run_dir.exists():
        identity_values, identity_issues = collect_identity_strings(run_dir, slot_map)
        issues.extend(identity_issues)
    else:
        issues.append("Identity evaluation skipped: run_dir missing or does not exist.")

    identity_condition: Optional[bool] = None
    identity_condition_issues: List[str] = []
    if identity_values:
        identity_condition, identity_condition_issues = evaluate_identity_condition(identity_values)
        issues.extend(identity_condition_issues)

    excluded = (genus_condition is False) and (identity_condition is False)
    filter_details = {
        "genus_condition_met": genus_condition,
        "identity_condition_met": identity_condition,
        "excluded": excluded,
    }
    return filter_details, issues


def compute_observed_vector(df: pd.DataFrame) -> np.ndarray:
    required_columns = {"mutType", "mutNum", "totalRootNum"}
    if missing := required_columns - set(df.columns):
        raise ValueError(f"TSV file is missing required columns: {', '.join(sorted(missing))}")

    df = df.copy()
    totals = df["totalRootNum"].to_numpy(dtype=float)
    mut_counts = df["mutNum"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        observed = np.divide(mut_counts, totals, out=np.zeros_like(mut_counts, dtype=float), where=totals != 0)
    return observed


def compute_logratio_vector(df: pd.DataFrame, observed: Optional[np.ndarray] = None) -> np.ndarray:
    if observed is None:
        observed = compute_observed_vector(df)
    mean_observed = np.mean(observed) if observed.size else 0.0
    adjusted_mean = max(mean_observed, EPS)
    ratios = (observed + EPS) / adjusted_mean
    with np.errstate(divide="ignore"):
        logratio = np.log2(ratios)
    return logratio


def load_tsv_vectors(
    tsv_path: Path, reference_order: Optional[List[str]]
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    df = pd.read_csv(tsv_path, sep="\t")
    if reference_order is None:
        # Sort by mutType to ensure consistent ordering of mutation types,
        # which is crucial for reproducible PCA/projection and for aligning
        # feature vectors across multiple samples.
        df_sorted = df.sort_values("mutType").reset_index(drop=True)
        order = df_sorted["mutType"].tolist()
        observed = compute_observed_vector(df_sorted)
        logratio = compute_logratio_vector(df_sorted, observed)
        return logratio, observed, order

    df_indexed = df.set_index("mutType")
    try:
        df_aligned = df_indexed.loc[reference_order].reset_index()
    except KeyError as exc:
        raise ValueError(f"TSV file {tsv_path.name} does not contain expected mutType values: {exc}") from exc
    observed = compute_observed_vector(df_aligned)
    logratio = compute_logratio_vector(df_aligned, observed)
    return logratio, observed, reference_order


def determine_clusters(data: np.ndarray, max_clusters: int, random_state: int) -> Tuple[int, List[float], List[float]]:
    sample_count = data.shape[0]
    if sample_count < 2:
        raise ValueError("At least two samples are required for clustering.")

    if sample_count == 2:
        # With only two samples the optimal clustering is trivially k=2, and silhouette
        # scores are undefined (n_labels must be <= n_samples - 1). Skip metric curves.
        return 2, [], []

    adjusted_max = min(max_clusters, sample_count - 1)
    if adjusted_max < 2:
        adjusted_max = 2

    inertia_values: List[float] = []
    silhouette_scores: List[float] = []

    for k in range(2, adjusted_max + 1):
        model = KMeans(n_clusters=k, random_state=random_state, n_init=10)
        labels = model.fit_predict(data)
        inertia_values.append(model.inertia_)
        silhouette_scores.append(silhouette_score(data, labels))

    best_k = silhouette_scores.index(max(silhouette_scores)) + 2
    return best_k, inertia_values, silhouette_scores


def plot_metrics(output_path: Path, inertia: List[float], silhouette: List[float]) -> None:
    if not inertia and not silhouette:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    if inertia:
        ks_inertia = list(range(2, 2 + len(inertia)))
        axes[0].plot(ks_inertia, inertia, "bo-")
        axes[0].set_xlabel("Number of Clusters")
        axes[0].set_ylabel("Inertia")
        axes[0].set_title("Elbow Method")
        axes[0].grid(True)
    else:
        axes[0].axis("off")
        axes[0].text(0.5, 0.5, "Insufficient samples\nfor inertia plot", ha="center", va="center")

    if silhouette:
        ks_silhouette = list(range(2, 2 + len(silhouette)))
        axes[1].plot(ks_silhouette, silhouette, "ro-")
        axes[1].set_xlabel("Number of Clusters")
        axes[1].set_ylabel("Silhouette Score")
        axes[1].set_title("Silhouette Scores")
        axes[1].grid(True)
    else:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "Silhouette undefined\n(n_samples too small)", ha="center", va="center")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_color_map(genus_labels: List[str]) -> Dict[str, Tuple[float, float, float, float]]:
    unique_genus = sorted(OrderedDict.fromkeys(genus_labels))
    cmap = plt.get_cmap("tab20")
    color_map: Dict[str, Tuple[float, float, float, float]] = {}
    for idx, genus in enumerate(unique_genus):
        color_map[genus] = cmap(idx % cmap.N)
    return color_map


def plot_pca_scatter(
    output_path: Path,
    pca_data: np.ndarray,
    species_labels: List[str],
    color_labels: List[str],
    cluster_labels: np.ndarray,
    label_name: str,
) -> None:
    color_map = build_color_map(color_labels)
    colors = [color_map[label] for label in color_labels]

    fig, ax = plt.subplots(figsize=(14, 12))
    ax.scatter(pca_data[:, 0], pca_data[:, 1], c=colors, alpha=0.8, s=80, edgecolors="k")

    texts: List[plt.Text] = []
    for idx, raw_text in enumerate(species_labels):
        clean_text = raw_text.replace("_", " ")
        text = ax.text(
            pca_data[idx, 0],
            pca_data[idx, 1],
            clean_text,
            fontsize=9,
            fontstyle="italic",
            ha="left",
            va="bottom",
        )
        texts.append(text)

    adjust_text(
        texts,
        x=pca_data[:, 0],
        y=pca_data[:, 1],
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.6, alpha=0.6, shrinkA=6, shrinkB=6),
        autoalign="xy",
        only_move={"points": "y", "text": "xy"},
        expand_text=(1.2, 1.6),
        expand_points=(1.2, 1.6),
        force_text=(0.6, 0.6),
        force_points=(0.4, 0.4),
        lim=1000,
    )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"K-means Clustering (k={len(set(cluster_labels))})")
    ax.grid(True)

    handles = []
    legend_title = label_name.capitalize() if label_name else "Classification"
    for label, color in color_map.items():
        handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color,
                markersize=10,
                label=label,
            )
        )
    ax.legend(handles=handles, title=legend_title, loc="best")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_placeholder_plot(output_path: Path, message: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
    ax.set_title("Filtered PCA")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_empty_results_csv(output_csv_path: Path, classification_levels: List[str]) -> None:
    columns = [
        "file",
        "accession",
        "short_name",
        "date",
        "metadata_genus",
        "species",
        "dataset",
        "cluster",
        "PC1",
        "PC2",
    ]
    for level in classification_levels:
        columns.append(f"classification_{level}")
    df = pd.DataFrame(columns=columns)
    df.to_csv(output_csv_path, index=False)


def save_results(
    output_csv_path: Path,
    files: List[Path],
    pca_data: np.ndarray,
    cluster_labels: np.ndarray,
    metadata_records: List[dict],
    classification_levels: List[str],
) -> None:
    rows: List[dict] = []
    for idx, file_path in enumerate(files):
        record = {
            "file": file_path.name,
            "accession": metadata_records[idx]["accession"],
            "short_name": metadata_records[idx]["short_name"],
            "date": metadata_records[idx]["date"],
            "metadata_genus": metadata_records[idx]["genus"],
            "species": metadata_records[idx]["species"],
            "dataset": metadata_records[idx]["dataset"],
            "cluster": int(cluster_labels[idx]),
            "PC1": float(pca_data[idx, 0]),
            "PC2": float(pca_data[idx, 1]),
        }
        for level in classification_levels:
            record[f"classification_{level}"] = metadata_records[idx]["classifications"].get(level, "Unknown")
        rows.append(record)

    df = pd.DataFrame(rows)
    df.to_csv(output_csv_path, index=False)


def main() -> None:
    args = parse_args()

    input_root = Path(args.input_dir).resolve()
    if not input_root.is_dir():
        print(f"Error: input directory does not exist: {input_root}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output_dir).resolve() if args.output_dir else input_root / "pca"
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_cache: Dict[Path, dict] = {}
    taxonomy_cache: Dict[Path, dict] = {}

    logratio_vectors: List[np.ndarray] = []
    observed_vectors: List[np.ndarray] = []
    metadata_records: List[dict] = []
    species_labels: List[str] = []
    source_files: List[Path] = []
    reference_order: Optional[List[str]] = None
    observed_levels: Set[str] = set()
    dataset_status: List[Tuple[str, Dict[str, Optional[bool]]]] = []
    dataset_sample_counts: Dict[str, int] = {}

    dataset_candidates = [
        entry for entry in sorted(input_root.iterdir()) if entry.is_dir() and not entry.name.isdigit()
    ]
    if not dataset_candidates:
        dataset_candidates = [input_root]

    for dataset_dir in dataset_candidates:
        dataset_name = dataset_dir.name
        latest_run = find_latest_run_dir(dataset_dir)
        if not latest_run:
            warn(f"{dataset_name}: No dated run directories found under {dataset_dir}.")
            continue

        manifest_path = latest_run / "metadata" / "metadata_manifest.json"
        if not manifest_path.exists():
            warn(f"{dataset_name}: metadata manifest not found at {manifest_path}.")
            continue

        try:
            _, manifest_lookup, _, slot_map = load_manifest_lookup(manifest_path)
        except Exception as exc:
            warn(f"{dataset_name}: Failed to load manifest {manifest_path}: {exc}")
            continue

        filter_info, filter_issues = evaluate_dataset_filter(slot_map, latest_run)
        for issue in filter_issues:
            warn(f"{dataset_name}: {issue}")
        dataset_status.append((dataset_name, filter_info))
        dataset_filter_flag = bool(filter_info.get("excluded"))

        tsv_files = sorted(latest_run.glob(args.tsv_pattern))
        if not tsv_files:
            warn(f"{dataset_name}: No TSV files matching pattern '{args.tsv_pattern}' in {latest_run}.")
            continue

        samples_before = len(logratio_vectors)
    for tsv_file in tsv_files:
            try:
        accession, short_name, date = parse_filename_parts(tsv_file)
            except ValueError as exc:
                warn(f"{dataset_name}: {exc}")
                continue

        entry_key = (accession, short_name)
        manifest_entry = manifest_lookup.get(entry_key)
        if manifest_entry is None:
                warn(
                    f"{dataset_name}: No manifest entry for accession {accession} / short_name {short_name} "
                    f"(file {tsv_file.name})."
            )
                continue

        metadata_path = manifest_entry["metadata_path"]
        taxonomy_path = manifest_entry["taxonomy_path"]

            try:
        metadata_json = load_json_cached(metadata_path, metadata_cache)
            except Exception as exc:
                warn(f"{dataset_name}: Failed to load metadata JSON {metadata_path}: {exc}")
                continue

            ensure_taxonomy_json(taxonomy_path, metadata_json, short_name, accession)
            try:
        taxonomy_json = load_json_cached(taxonomy_path, taxonomy_cache)
            except Exception as exc:
                warn(f"{dataset_name}: Failed to load taxonomy JSON {taxonomy_path}: {exc}")
                taxonomy_json = {}

        genus, species_label = extract_genus(metadata_json, short_name, accession)
        classification_map = extract_classification_map(taxonomy_json, f"{short_name} ({accession})")
        observed_levels.update(classification_map.keys())

            try:
                logratio_vec, observed_vec, reference_order = load_tsv_vectors(tsv_file, reference_order)
            except Exception as exc:
                warn(f"{dataset_name}: Failed to parse TSV {tsv_file.name}: {exc}")
                continue

            logratio_vectors.append(logratio_vec)
            observed_vectors.append(observed_vec)
        species_labels.append(species_label)
            source_files.append(tsv_file)
        metadata_records.append(
            {
                "accession": accession,
                "short_name": short_name,
                "date": date,
                "genus": genus,
                "species": species_label,
                "classifications": classification_map,
                    "dataset": dataset_name,
                    "filter_excluded": dataset_filter_flag,
                }
            )

        added = len(logratio_vectors) - samples_before
        if added:
            dataset_sample_counts[dataset_name] = added
        else:
            warn(f"{dataset_name}: No usable TSV files after validation; dataset skipped.")

    if len(logratio_vectors) < 2:
        print("Error: At least two valid TSV samples are required for PCA.", file=sys.stderr)
        sys.exit(1)

    ordered_levels = [level for level in CLASSIFICATION_ORDER if level in observed_levels]
    remaining_levels = sorted(observed_levels - set(ordered_levels))
    ordered_levels.extend(remaining_levels)
    if not ordered_levels:
        ordered_levels = ["taxonomy"]

    plot_levels = [level for level in ordered_levels if level != "species"]
    if not plot_levels:
        print("Warning: Only species-level classifications detected; skipping classification plots.", file=sys.stderr)
    sample_count = len(logratio_vectors)
    full_indices = list(range(sample_count))

    def run_pca_for_indices(
        vectors: List[np.ndarray],
        indices: List[int],
        metrics_path: Path,
        csv_path: Path,
        plot_prefix: str,
    ) -> Tuple[int, int]:
        subset_vectors = [vectors[i] for i in indices]
        data_matrix = np.vstack(subset_vectors)
        scaler = StandardScaler()
        scaled_subset = scaler.fit_transform(data_matrix)

        pca = PCA(n_components=2)
        pca_data = pca.fit_transform(scaled_subset)

        best_k_value, inertia_values, silhouette_scores = determine_clusters(pca_data, args.max_clusters, args.random_state)
        plot_metrics(metrics_path, inertia_values, silhouette_scores)

        final_model = KMeans(n_clusters=best_k_value, random_state=args.random_state, n_init=10)
        cluster_labels = final_model.fit_predict(pca_data)

        subset_species = [species_labels[i] for i in indices]
        subset_metadata = [metadata_records[i] for i in indices]
        subset_files = [source_files[i] for i in indices]

        if plot_levels:
            for level in plot_levels:
                labels = [record["classifications"].get(level, "Unknown") for record in subset_metadata]
                sanitized_level = level.replace(" ", "_")
                output_path = output_dir / f"{plot_prefix}_{sanitized_level}.png"
                plot_pca_scatter(output_path, pca_data, subset_species, labels, cluster_labels, level)

        save_results(csv_path, subset_files, pca_data, cluster_labels, subset_metadata, ordered_levels)
        return best_k_value, len(subset_files)

    filtered_indices = [idx for idx, record in enumerate(metadata_records) if not record.get("filter_excluded")]

    def process_variant(vectors: List[np.ndarray], variant_label: str) -> Dict[str, object]:
        metrics_path = output_dir / f"clustering_metrics_{variant_label}.png"
        csv_path = output_dir / f"clustering_results_{variant_label}.csv"
        plot_prefix = f"pca_2d_clusters_{variant_label}"
        best_k_value, processed_count = run_pca_for_indices(
            vectors, full_indices, metrics_path, csv_path, plot_prefix
        )

        result: Dict[str, object] = {
            "label": variant_label,
            "best_k": best_k_value,
            "processed": processed_count,
            "filtered_best_k": None,
            "filtered_processed": 0,
            "filtered_message": None,
        }

        if len(filtered_indices) >= 2:
            filtered_best_k, filtered_processed = run_pca_for_indices(
                vectors,
                filtered_indices,
                output_dir / f"clustering_metrics_{variant_label}_filtered.png",
                output_dir / f"clustering_results_{variant_label}_filtered.csv",
                f"pca_2d_clusters_{variant_label}_filtered",
            )
            result["filtered_best_k"] = filtered_best_k
            result["filtered_processed"] = filtered_processed
        else:
            if filtered_indices:
                filtered_message = (
                    f"{variant_label} filtered PCA requires at least two samples; skipping filtered computation."
                )
    else:
                filtered_message = (
                    f"{variant_label} filtered PCA skipped: no samples remained after applying the filter."
                )
            write_placeholder_plot(
                output_dir / f"clustering_metrics_{variant_label}_filtered.png", filtered_message
            )
            write_empty_results_csv(output_dir / f"clustering_results_{variant_label}_filtered.csv", ordered_levels)
            if plot_levels:
        for level in plot_levels:
            sanitized_level = level.replace(" ", "_")
                    placeholder_path = (
                        output_dir / f"pca_2d_clusters_{variant_label}_filtered_{sanitized_level}.png"
                    )
                    write_placeholder_plot(placeholder_path, filtered_message)
            result["filtered_message"] = filtered_message

        return result

    variant_results = [
        process_variant(logratio_vectors, "logratio"),
        process_variant(observed_vectors, "observed"),
    ]

    print("PCA clustering completed successfully.")
    print(f"Datasets processed: {len(dataset_sample_counts)}")
    print(f"Total samples processed: {sample_count}")
    for entry in variant_results:
        label = entry["label"]
        print(
            f"{label.capitalize()} variant — processed: {entry['processed']}, optimal k: {entry['best_k']}"
        )
        filtered_processed = entry["filtered_processed"]
        filtered_best_k = entry["filtered_best_k"]
        filtered_message = entry["filtered_message"]
        if filtered_processed:
            print(f"  Filtered processed files: {filtered_processed}, optimal k: {filtered_best_k}")
        elif filtered_message:
            print(f"  {filtered_message}")
    if plot_levels:
        print("Generated classification plots for levels:", ", ".join(plot_levels))
    else:
        print("No classification plots generated (species-only classifications).")
    if dataset_status:
        print("Filter evaluation summary:")
        for dataset_name, filter_info in dataset_status:
            print(
                f"  {dataset_name}: excluded={filter_info.get('excluded')} "
                f"genus_condition={filter_info.get('genus_condition_met')} "
                f"identity_condition={filter_info.get('identity_condition_met')}"
            )
    print(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

