#!/usr/bin/env python3
"""
Authors: Simon Kraft, Joshua Payne, El Sall
CPSC 482 - Data Structures II

Build report tables from benchmark and optimization CSV files.

Outputs:
  out/tables/problem_stats.md
  out/tables/optimization_stage_comparison.md
  out/tables/algorithm_runtime_summary.md
"""

import os
from pathlib import Path
import pandas as pd


def write_md_table(df: pd.DataFrame, path: Path, title: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    headers = [str(c) for c in df.columns.tolist()]
    rows = [[str(v) for v in row] for row in df.values.tolist()]
    with path.open("w", encoding="utf-8") as f:
        f.write(f"# {title}\n\n")
        f.write("| " + " | ".join(headers) + " |\n")
        f.write("|" + "|".join(["---"] * len(headers)) + "|\n")
        for row in rows:
            f.write("| " + " | ".join(row) + " |\n")
        f.write("\n")


def build_problem_stats(opt_df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "matrix", "n", "nnz", "density_pct",
        "fbm_base", "tfbd_base"
    ]
    out = opt_df[cols].copy().sort_values("n")
    out = out.rename(columns={
        "density_pct": "density(%)",
        "fbm_base": "FBM(base)",
        "tfbd_base": "TFBD(base)",
    })
    out["density(%)"] = out["density(%)"].map(lambda x: f"{x:.4f}")
    return out


def build_optimization_table(opt_df: pd.DataFrame) -> pd.DataFrame:
    sorted_df = opt_df.sort_values("n" if "n" in opt_df.columns else "matrix")
    cols = [
        "matrix",
        "fbm_base", "fbm_topo", "fbm_tearing", "fbm_banding",
        "tfbd_base", "tfbd_topo", "tfbd_tearing", "tfbd_banding",
        "fbm_gain_banding_pct", "tfbd_gain_banding_pct",
    ]
    out = sorted_df[cols].copy()
    out = out.rename(columns={
        "fbm_gain_banding_pct": "FBM gain final(%)",
        "tfbd_gain_banding_pct": "TFBD gain final(%)",
    })
    out["FBM gain final(%)"] = out["FBM gain final(%)"].map(lambda x: f"{x:.2f}")
    out["TFBD gain final(%)"] = out["TFBD gain final(%)"].map(lambda x: f"{x:.2f}")
    return out


def build_runtime_summary(bench_df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        bench_df.groupby("algorithm", as_index=False)
        .agg(
            mean_us=("mean_us", "mean"),
            median_us=("median_us", "mean"),
            stddev_us=("stddev_us", "mean"),
            mean_ns_per_vertex=("ns_per_vertex", "mean"),
            mean_ns_per_edge=("ns_per_edge", "mean"),
        )
        .sort_values("mean_us")
    )
    for col in ["mean_us", "median_us", "stddev_us", "mean_ns_per_vertex", "mean_ns_per_edge"]:
        summary[col] = summary[col].map(lambda x: f"{x:.2f}")
    return summary


def main():
    root = Path(__file__).resolve().parent.parent
    benchmark_csv = root / "out" / "benchmark.csv"
    optimization_csv = root / "out" / "optimization_metrics.csv"

    if not benchmark_csv.exists() or not optimization_csv.exists():
        print("Missing CSV input files.")
        print("Run ./out/main --bench first.")
        raise SystemExit(1)

    bench_df = pd.read_csv(benchmark_csv)
    opt_df = pd.read_csv(optimization_csv)
    # Align runtime averages with report/plots (10-matrix set).
    bench_df = bench_df[bench_df["matrix"] != "easy-example"].copy()
    opt_df = opt_df[opt_df["matrix"] != "easy-example"].copy()

    tables_dir = root / "out" / "tables"

    write_md_table(
        build_problem_stats(opt_df),
        tables_dir / "problem_stats.md",
        "Problem Statistics",
    )
    write_md_table(
        build_optimization_table(opt_df),
        tables_dir / "optimization_stage_comparison.md",
        "Optimization Stage Comparison",
    )
    write_md_table(
        build_runtime_summary(bench_df),
        tables_dir / "algorithm_runtime_summary.md",
        "Algorithm Runtime Summary",
    )

    print(f"Tables written to: {tables_dir}")
    for p in sorted(tables_dir.glob("*.md")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
