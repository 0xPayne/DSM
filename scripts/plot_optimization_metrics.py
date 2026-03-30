#!/usr/bin/env python3
"""
Generate report-ready plots from out/optimization_metrics.csv.

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


def grouped_reduction_plot(df: pd.DataFrame, metric_prefix: str, out_path: str, title: str, ylabel: str):
    matrices = df.sort_values("n")["matrix"].tolist()
    x = np.arange(len(matrices))
    width = 0.25

    topo = df.set_index("matrix")[f"{metric_prefix}_gain_topo_pct"].reindex(matrices).values
    tear = df.set_index("matrix")[f"{metric_prefix}_gain_tearing_pct"].reindex(matrices).values
    band = df.set_index("matrix")[f"{metric_prefix}_gain_banding_pct"].reindex(matrices).values

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
    summary = {
        "FBM reduction (%)": [
            df["fbm_gain_topo_pct"].mean(),
            df["fbm_gain_tearing_pct"].mean(),
            df["fbm_gain_banding_pct"].mean(),
        ],
        "TFBD reduction (%)": [
            df["tfbd_gain_topo_pct"].mean(),
            df["tfbd_gain_tearing_pct"].mean(),
            df["tfbd_gain_banding_pct"].mean(),
        ],
    }

    stages = ["TopoOnly", "Tearing", "Banding"]
    x = np.arange(len(stages))
    width = 0.35

    fig, ax = plt.subplots()
    ax.bar(x - width / 2, summary["FBM reduction (%)"], width, label="FBM", color="#0891b2")
    ax.bar(x + width / 2, summary["TFBD reduction (%)"], width, label="TFBD", color="#ca8a04")
    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    ax.set_ylabel("Average reduction vs baseline (%)")
    ax.set_title("Average Optimization Gain by Stage")
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
        title="FBM Reduction by Reordering Stage",
        ylabel="Reduction vs baseline (%)",
    )
    grouped_reduction_plot(
        df=df,
        metric_prefix="tfbd",
        out_path=os.path.join(out_dir, "optimization_tfbd_reduction.png"),
        title="TFBD Reduction by Reordering Stage",
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
