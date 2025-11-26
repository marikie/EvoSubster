#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

rcParams["font.family"] = "MS Gothic"

EPS = 1e-12
DEFAULT_TSV_PATTERN = "*.tsv"


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run 2D PCA clustering on logRatio vectors derived from substitution TSV files."
    )
    parser.add_argument(
        "--tsv-dir",
        required=True,
        help=(
            "Directory containing substitution TSV files (e.g. outputs of generate_tsv_files.sh). "
            f"Files should follow the pattern '{DEFAULT_TSV_PATTERN}' unless --tsv-pattern is provided."
        ),
    )
    parser.add_argument(
        "--output-dir",
        help="Directory where PCA outputs will be written. Defaults to <tsv-dir>/pca_logratio_output.",
    )
    parser.add_argument(
        "--tsv-pattern",
        default=DEFAULT_TSV_PATTERN,
        help=f"Glob pattern to select TSV files within the TSV directory (default: {DEFAULT_TSV_PATTERN}).",
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
            "\nExample:\n  python scripts/analysis/pca_logratio_2d.py "
            "--tsv-dir /path/to/tsv_dir --output-dir /path/to/output_dir",
            file=sys.stderr,
        )
        raise SystemExit(1)
    return parser.parse_args(argv)


def collect_tsv_files(tsv_dir: Path, pattern: str) -> List[Path]:
    files = sorted(tsv_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No TSV files matching pattern '{pattern}' were found in {tsv_dir}")
    return files


def parse_filename_parts(tsv_path: Path) -> Tuple[str, str, str]:
    name_parts = tsv_path.stem.split("_")
    if len(name_parts) < 3:
        raise ValueError(f"Unexpected filename format for TSV file: {tsv_path.name}")
    accession = name_parts[0]
    short_name = name_parts[1]
    date = name_parts[2]
    return accession, short_name, date


def run_datasets_cli(accession: str) -> dict:
    try:
        proc = subprocess.run(
            ["datasets", "summary", "genome", "accession", accession],
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        raise RuntimeError(f"Failed to fetch metadata for accession {accession}: {exc}") from exc

    output = proc.stdout.strip()
    if not output:
        raise RuntimeError(f"Empty metadata response for accession {accession}")

    try:
        return json.loads(output)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Invalid JSON received for accession {accession}: {exc}") from exc


def ensure_metadata(tsv_path: Path, accession: str) -> dict:
    metadata_path = tsv_path.with_name(f"{accession}.json")
    if metadata_path.exists():
        with metadata_path.open("r", encoding="utf-8") as handle:
            return json.load(handle)

    metadata = run_datasets_cli(accession)
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)
    return metadata


def extract_genus(metadata: dict, short_name: str, accession: str) -> Tuple[str, str]:
    reports = metadata.get("reports", [])
    if reports:
        organism = reports[0].get("organism", {})
        organism_name = organism.get("organism_name")
        if organism_name:
            genus = organism_name.split()[0]
            label = f"{organism_name.replace(' ', '_')} ({short_name})"
            return genus, label

    fallback_label = f"{short_name} ({accession})"
    return "Unknown", fallback_label


def compute_logratio_vector(df: pd.DataFrame) -> np.ndarray:
    required_columns = {"mutType", "mutNum", "totalRootNum"}
    if missing := required_columns - set(df.columns):
        raise ValueError(f"TSV file is missing required columns: {', '.join(sorted(missing))}")

    df = df.copy()
    totals = df["totalRootNum"].to_numpy(dtype=float)
    mut_counts = df["mutNum"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        observed = np.divide(mut_counts, totals, out=np.zeros_like(mut_counts, dtype=float), where=totals != 0)
    mean_observed = np.mean(observed) if observed.size else 0.0
    # Use max(mean_observed, EPS) to avoid division by zero or producing extremely large ratios
    # when mean_observed is zero or extremely close to zero. Ensures numerical stability.
    adjusted_mean = max(mean_observed, EPS)
    ratios = (observed + EPS) / adjusted_mean
    with np.errstate(divide="ignore"):
        logratio = np.log2(ratios)
    return logratio


def load_tsv_logratio(tsv_path: Path, reference_order: Optional[List[str]]) -> Tuple[np.ndarray, List[str]]:
    df = pd.read_csv(tsv_path, sep="\t")
    if reference_order is None:
        # Sort by mutType to ensure consistent ordering of mutation types,
        # which is crucial for reproducible PCA/projection and for aligning
        # feature vectors across multiple samples.
        df_sorted = df.sort_values("mutType").reset_index(drop=True)
        order = df_sorted["mutType"].tolist()
        vector = compute_logratio_vector(df_sorted)
        return vector, order

    df_indexed = df.set_index("mutType")
    try:
        df_aligned = df_indexed.loc[reference_order].reset_index()
    except KeyError as exc:
        raise ValueError(f"TSV file {tsv_path.name} does not contain expected mutType values: {exc}") from exc
    vector = compute_logratio_vector(df_aligned)
    return vector, reference_order


def determine_clusters(data: np.ndarray, max_clusters: int, random_state: int) -> Tuple[int, List[float], List[float]]:
    sample_count = data.shape[0]
    if sample_count < 2:
        raise ValueError("At least two samples are required for clustering.")

    adjusted_max = min(max_clusters, sample_count)
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


def plot_metrics(output_dir: Path, inertia: List[float], silhouette: List[float]) -> None:
    ks = list(range(2, 2 + len(inertia)))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(ks, inertia, "bo-")
    axes[0].set_xlabel("Number of Clusters")
    axes[0].set_ylabel("Inertia")
    axes[0].set_title("Elbow Method")
    axes[0].grid(True)

    axes[1].plot(ks, silhouette, "ro-")
    axes[1].set_xlabel("Number of Clusters")
    axes[1].set_ylabel("Silhouette Score")
    axes[1].set_title("Silhouette Scores")
    axes[1].grid(True)

    plt.tight_layout()
    fig.savefig(output_dir / "clustering_metrics.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_color_map(genus_labels: List[str]) -> Dict[str, Tuple[float, float, float, float]]:
    unique_genus = sorted(OrderedDict.fromkeys(genus_labels))
    cmap = plt.get_cmap("tab20")
    color_map: Dict[str, Tuple[float, float, float, float]] = {}
    for idx, genus in enumerate(unique_genus):
        color_map[genus] = cmap(idx % cmap.N)
    return color_map


def plot_pca_scatter(
    output_dir: Path,
    pca_data: np.ndarray,
    labels: List[str],
    genus_labels: List[str],
    cluster_labels: np.ndarray,
) -> None:
    color_map = build_color_map(genus_labels)
    colors = [color_map[genus] for genus in genus_labels]

    fig, ax = plt.subplots(figsize=(14, 12))
    ax.scatter(pca_data[:, 0], pca_data[:, 1], c=colors, alpha=0.8, s=80, edgecolors="k")

    for idx, text in enumerate(labels):
        ax.annotate(text, (pca_data[idx, 0], pca_data[idx, 1]), xytext=(5, 5), textcoords="offset points", fontsize=8)

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"K-means Clustering (k={len(set(cluster_labels))})")
    ax.grid(True)

    handles = []
    for genus, color in color_map.items():
        handles.append(plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=10, label=genus))
    ax.legend(handles=handles, title="Genus", loc="best")

    fig.savefig(output_dir / "pca_2d_clusters.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_results(
    output_dir: Path,
    files: List[Path],
    pca_data: np.ndarray,
    cluster_labels: np.ndarray,
    metadata_records: List[dict],
) -> None:
    rows: List[dict] = []
    for idx, file_path in enumerate(files):
        record = {
            "file": file_path.name,
            "accession": metadata_records[idx]["accession"],
            "short_name": metadata_records[idx]["short_name"],
            "date": metadata_records[idx]["date"],
            "genus": metadata_records[idx]["genus"],
            "cluster": int(cluster_labels[idx]),
            "PC1": float(pca_data[idx, 0]),
            "PC2": float(pca_data[idx, 1]),
        }
        rows.append(record)

    df = pd.DataFrame(rows)
    df.to_csv(output_dir / "clustering_results.csv", index=False)


def main() -> None:
    args = parse_args()

    tsv_dir = Path(args.tsv_dir).resolve()
    if not tsv_dir.is_dir():
        print(f"Error: TSV directory does not exist: {tsv_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output_dir).resolve() if args.output_dir else tsv_dir / "pca_logratio_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    tsv_files = collect_tsv_files(tsv_dir, args.tsv_pattern)

    vectors: List[np.ndarray] = []
    metadata_records: List[dict] = []
    labels: List[str] = []
    genus_labels: List[str] = []
    reference_order: Optional[List[str]] = None

    for tsv_file in tsv_files:
        accession, short_name, date = parse_filename_parts(tsv_file)
        metadata_json = ensure_metadata(tsv_file, accession)
        genus, label = extract_genus(metadata_json, short_name, accession)

        vector, reference_order = load_tsv_logratio(tsv_file, reference_order)

        vectors.append(vector)
        labels.append(label)
        genus_labels.append(genus)
        metadata_records.append(
            {
                "accession": accession,
                "short_name": short_name,
                "date": date,
                "genus": genus,
            }
        )

    data_matrix = np.vstack(vectors)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data_matrix)

    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(scaled_data)

    best_k, inertia_values, silhouette_scores = determine_clusters(pca_data, args.max_clusters, args.random_state)
    plot_metrics(output_dir, inertia_values, silhouette_scores)

    final_model = KMeans(n_clusters=best_k, random_state=args.random_state, n_init=10)
    cluster_labels = final_model.fit_predict(pca_data)

    plot_pca_scatter(output_dir, pca_data, labels, genus_labels, cluster_labels)
    save_results(output_dir, tsv_files, pca_data, cluster_labels, metadata_records)

    print("PCA clustering completed successfully.")
    print(f"Processed files: {len(tsv_files)}")
    print(f"Optimal clusters (k): {best_k}")
    print(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

