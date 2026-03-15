#!/usr/bin/env python3
"""
Generate publication-quality plots from DSM benchmark CSV output.

Usage:
    python3 scripts/plot_benchmarks.py [path/to/benchmark.csv]

Default CSV path: out/benchmark.csv
Output:  out/plots/*.png
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ---------------------------------------------------------------------------
# Style: clean, publication-ready
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "figure.figsize":    (10, 6),
    "figure.dpi":        150,
    "font.family":       "serif",
    "font.size":         11,
    "axes.titlesize":    13,
    "axes.labelsize":    12,
    "legend.fontsize":   10,
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "axes.grid":         True,
    "grid.alpha":        0.3,
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

COLORS = {
    "Tarjan":       "#2563eb",
    "Kosaraju":     "#dc2626",
    "Condensation": "#16a34a",
    "TopoSort":     "#9333ea",
    "Permute":      "#ea580c",
    "FBM-Count":    "#0891b2",
    "TFBD-Sum":     "#ca8a04",
}

def load_data(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["ve_sum"] = df["n"] + df["nnz"]  # V + E for complexity analysis
    return df


def plot_scaling(df: pd.DataFrame, out_dir: str):
    """Log-log: mean runtime vs (V + E) per algorithm.
    For O(V+E) algorithms, expect slope ~1."""
    fig, ax = plt.subplots()

    for algo, group in df.groupby("algorithm"):
        g = group.sort_values("ve_sum")
        color = COLORS.get(algo, "#6b7280")
        ax.loglog(g["ve_sum"], g["mean_us"], "o-", label=algo,
                  color=color, markersize=5, linewidth=1.5)

    # Reference line: slope = 1 (linear)
    xr = np.array([df["ve_sum"].min(), df["ve_sum"].max()])
    scale = df.loc[df["ve_sum"].idxmax(), "mean_us"] / df["ve_sum"].max()
    ax.loglog(xr, xr * scale, "--", color="#9ca3af", linewidth=1, label="O(V+E) ref")

    ax.set_xlabel("V + E  (vertices + edges)")
    ax.set_ylabel("Mean runtime (μs)")
    ax.set_title("Algorithm Scaling: Runtime vs Problem Size")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "scaling.png"))
    plt.close(fig)


def plot_throughput(df: pd.DataFrame, out_dir: str):
    """ns per (V+E) vs problem size. Flat = linear scaling confirmed."""
    fig, ax = plt.subplots()

    for algo, group in df.groupby("algorithm"):
        g = group.sort_values("ve_sum")
        throughput = g["mean_us"].values * 1000.0 / g["ve_sum"].values
        color = COLORS.get(algo, "#6b7280")
        ax.semilogx(g["ve_sum"], throughput, "s-", label=algo,
                     color=color, markersize=5, linewidth=1.5)

    ax.set_xlabel("V + E  (vertices + edges)")
    ax.set_ylabel("ns / (V + E)")
    ax.set_title("Throughput: Amortized Cost per Graph Element")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "throughput.png"))
    plt.close(fig)


def plot_per_matrix_bars(df: pd.DataFrame, out_dir: str):
    """Grouped bar chart: runtime per matrix, grouped by algorithm."""
    matrices = df.sort_values("n")["matrix"].unique()
    algorithms = df["algorithm"].unique()
    x = np.arange(len(matrices))
    width = 0.8 / len(algorithms)

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, algo in enumerate(algorithms):
        sub = df[df["algorithm"] == algo].set_index("matrix")
        vals = [sub.loc[m, "mean_us"] if m in sub.index else 0 for m in matrices]
        color = COLORS.get(algo, "#6b7280")
        ax.bar(x + i * width, vals, width, label=algo, color=color, alpha=0.85)

    ax.set_xticks(x + width * (len(algorithms) - 1) / 2)
    ax.set_xticklabels(matrices, rotation=45, ha="right")
    ax.set_ylabel("Mean runtime (μs)")
    ax.set_title("Per-Matrix Algorithm Comparison")
    ax.set_yscale("log")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "per_matrix.png"))
    plt.close(fig)


def plot_variability(df: pd.DataFrame, out_dir: str):
    """Coefficient of variation (stddev/mean) per measurement.
    Low CV = stable, reliable timing."""
    df_cv = df.copy()
    df_cv["cv"] = df_cv["stddev_us"] / df_cv["mean_us"].replace(0, np.nan)
    df_cv = df_cv.dropna(subset=["cv"])

    fig, ax = plt.subplots()

    for algo, group in df_cv.groupby("algorithm"):
        g = group.sort_values("ve_sum")
        color = COLORS.get(algo, "#6b7280")
        ax.semilogx(g["ve_sum"], g["cv"] * 100, "o-", label=algo,
                     color=color, markersize=5, linewidth=1.5)

    ax.set_xlabel("V + E  (vertices + edges)")
    ax.set_ylabel("Coefficient of Variation (%)")
    ax.set_title("Timing Stability Across Problem Sizes")
    ax.axhline(y=5, color="#ef4444", linestyle=":", linewidth=1, alpha=0.5, label="5% threshold")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "variability.png"))
    plt.close(fig)


def plot_individual_runs(df: pd.DataFrame, out_dir: str):
    """Strip plot of all individual run times for each matrix (Tarjan only)."""
    tarjan = df[df["algorithm"] == "Tarjan"].sort_values("n")
    run_cols = [c for c in tarjan.columns if c.startswith("run") and c.endswith("_us")]
    if not run_cols:
        return

    fig, ax = plt.subplots(figsize=(12, 5))
    positions = []
    labels = []

    for i, (_, row) in enumerate(tarjan.iterrows()):
        runs = row[run_cols].values.astype(float)
        jitter = np.random.normal(0, 0.08, size=len(runs))
        ax.scatter([i + j for j in jitter], runs,
                   color=COLORS["Tarjan"], alpha=0.7, s=30, zorder=3)
        ax.scatter([i], [row["mean_us"]], color="#dc2626",
                   marker="D", s=50, zorder=4, edgecolors="white", linewidths=0.5)
        positions.append(i)
        labels.append(row["matrix"])

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Runtime (μs)")
    ax.set_yscale("log")
    ax.set_title("Tarjan SCC: Individual Run Distribution (◆ = mean)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "run_distribution.png"))
    plt.close(fig)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    default_csv = os.path.join(project_root, "out", "benchmark.csv")

    csv_path = sys.argv[1] if len(sys.argv) > 1 else default_csv

    if not os.path.exists(csv_path):
        print(f"CSV not found: {csv_path}")
        print("Run  ./out/main --bench  first to generate benchmark data.")
        sys.exit(1)

    out_dir = os.path.join(project_root, "out", "plots")
    os.makedirs(out_dir, exist_ok=True)

    df = load_data(csv_path)
    print(f"Loaded {len(df)} records from {csv_path}")
    print(f"Matrices: {df['matrix'].nunique()}, Algorithms: {df['algorithm'].nunique()}")

    plot_scaling(df, out_dir)
    plot_throughput(df, out_dir)
    plot_per_matrix_bars(df, out_dir)
    plot_variability(df, out_dir)
    plot_individual_runs(df, out_dir)

    print(f"\nPlots saved to: {out_dir}/")
    for f in sorted(os.listdir(out_dir)):
        if f.endswith(".png"):
            print(f"  {f}")


if __name__ == "__main__":
    main()
