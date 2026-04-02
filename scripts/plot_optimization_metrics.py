#!/usr/bin/env python3
"""
Authors: Simon Kraft, Joshua Payne, El Sall
CPSC 482 - Data Structures II

Generate report-ready plots from out/optimization_metrics.csv.

Per-matrix grouped bars use the CSV columns *_gain_{topo,tearing,banding}_pct:
each bar is cumulative reduction vs the baseline ordering (same definitions as
benchmark export and optimization_summary.csv).

Usage:
    python3 scripts/plot_optimization_metrics.py [path/to/optimization_metrics.csv]

Default CSV path: out/optimization_metrics.csv
Output: out/plots/optimization_*.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


plt.rcParams.update({
    "figure.figsize": (10, 6),
    "figure.dpi": 150,
    "font.family": "serif",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

STAGE_COLORS = {
    "TopoOnly": "#2563eb",
    "Tearing": "#ea580c",
    "Banding": "#16a34a",
}


def _cumulative_pct_cols(df: pd.DataFrame, kind: str):
    """Return (topo%, tear%, band%) vs baseline from precomputed CSV columns."""
    return (
        df[f"{kind}_gain_topo_pct"].to_numpy(),
        df[f"{kind}_gain_tearing_pct"].to_numpy(),
        df[f"{kind}_gain_banding_pct"].to_numpy(),
    )


def grouped_reduction_plot(df: pd.DataFrame, metric_prefix: str, out_path: str, title: str, ylabel: str):
    ordered = df.sort_values("n")
    matrices = ordered["matrix"].tolist()
    topo, tear, band = _cumulative_pct_cols(ordered, metric_prefix)

    x = np.arange(len(matrices))
    width = 0.25

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(x - width, topo, width, label="TopoOnly", color=STAGE_COLORS["TopoOnly"])
    ax.bar(x, tear, width, label="Tearing", color=STAGE_COLORS["Tearing"])
    ax.bar(x + width, band, width, label="Banding", color=STAGE_COLORS["Banding"])

    ax.set_xticks(x)
    ax.set_xticklabels(matrices, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def cumulative_summary_plot(df: pd.DataFrame, out_path: str):
    fbm_topo, fbm_tear, fbm_band = _cumulative_pct_cols(df, "fbm")
    tfbd_topo, tfbd_tear, tfbd_band = _cumulative_pct_cols(df, "tfbd")
    summary = {
        "FBM (%)": [fbm_topo.mean(), fbm_tear.mean(), fbm_band.mean()],
        "TFBD (%)": [tfbd_topo.mean(), tfbd_tear.mean(), tfbd_band.mean()],
    }

    stages = ["TopoOnly", "Tearing", "Banding"]
    x = np.arange(len(stages))
    width = 0.35

    fig, ax = plt.subplots()
    ax.bar(x - width / 2, summary["FBM (%)"], width, label="FBM", color="#0891b2")
    ax.bar(x + width / 2, summary["TFBD (%)"], width, label="TFBD", color="#ca8a04")
    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    ax.set_ylabel("Mean reduction vs baseline (%)")
    ax.set_title("Average reduction by stage (10 matrices)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def absolute_gain_scatter(df: pd.DataFrame, out_path: str):
    fig, ax = plt.subplots()
    ax.scatter(df["n"], df["fbm_gain_banding_abs"], color="#16a34a", label="FBM absolute gain", s=35)
    ax.scatter(df["n"], df["tfbd_gain_banding_abs"], color="#dc2626", label="TFBD absolute gain", s=35)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Matrix size n (log scale)")
    ax.set_ylabel("Absolute gain at Banding stage (log scale)")
    ax.set_title("Final Absolute Gain vs Matrix Size")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    default_csv = os.path.join(project_root, "out", "optimization_metrics.csv")
    csv_path = sys.argv[1] if len(sys.argv) > 1 else default_csv

    if not os.path.exists(csv_path):
        print(f"CSV not found: {csv_path}")
        print("Run ./out/main --bench first to generate optimization metrics.")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    # Keep optimization plots aligned with the 10-matrix benchmark set.
    df = df[df["matrix"] != "easy-example"].copy()
    out_dir = os.path.join(project_root, "out", "plots")
    os.makedirs(out_dir, exist_ok=True)

    grouped_reduction_plot(
        df=df,
        metric_prefix="fbm",
        out_path=os.path.join(out_dir, "optimization_fbm_reduction.png"),
        title="FBM reduction by reordering stage",
        ylabel="Reduction vs baseline (%)",
    )
    grouped_reduction_plot(
        df=df,
        metric_prefix="tfbd",
        out_path=os.path.join(out_dir, "optimization_tfbd_reduction.png"),
        title="TFBD reduction by reordering stage",
        ylabel="Reduction vs baseline (%)",
    )
    cumulative_summary_plot(
        df=df,
        out_path=os.path.join(out_dir, "optimization_stage_average.png"),
    )
    absolute_gain_scatter(
        df=df,
        out_path=os.path.join(out_dir, "optimization_absolute_gain_vs_size.png"),
    )

    print(f"Loaded {len(df)} matrices from {csv_path}")
    print(f"Plots saved to: {out_dir}/")
    for file_name in sorted(os.listdir(out_dir)):
        if file_name.startswith("optimization_") and file_name.endswith(".png"):
            print(f"  {file_name}")


if __name__ == "__main__":
    main()
