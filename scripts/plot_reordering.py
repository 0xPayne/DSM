#!/usr/bin/env python3
"""
Plot sparse matrices before and after permutation.

Modes:
  1. Single matrix (positional arg):
       python scripts/plot_reordering.py data/easy-example.mtx
       python scripts/plot_reordering.py data/easy-example.mtx --perm perm.npy --output fig.png
       python scripts/plot_reordering.py data/easy-example.mtx --perm perm.txt --perm-format old2new

  2. Multi-matrix comparison (outputs one PNG per matrix):
       python scripts/plot_reordering.py --compare California GD99_c Tina_AskCal wb-cs-stanford \
           --markersizes 0.25 1.0 5.0 0.25 --output-dir report/figures/

Permutation formats supported (single matrix mode):
- "new2old" (default): array p where p[new_index] = old_index (0-based)
- "old2new": array q where q[old_index] = new_index; converted via np.argsort

SCC block boxes are drawn on permuted plots if out/sccs/<stem>_sccs.txt exists.
"""

import argparse
import numpy as np
import scipy.sparse as sp
from scipy.io import mmread
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path


# ── helpers ───────────────────────────────────────────────────────────────────

def load_permutation(path, n, fmt="new2old", one_based=False):
    perm = np.load(path) if path.endswith(".npy") else np.loadtxt(path, dtype=int)
    perm = perm.astype(int)
    if one_based:
        perm -= 1
    if perm.size != n:
        raise ValueError(f"Permutation length {perm.size} != matrix size {n}")
    if fmt == "old2new":
        perm = np.argsort(perm)
    elif fmt != "new2old":
        raise ValueError("Unknown perm format: choose 'new2old' or 'old2new'")
    return perm


def load_scc_blocks(sccs_path):
    with open(sccs_path) as f:
        return [int(l.strip()) for l in f if l.strip()]


def auto_markersize(n):
    if n < 100:   return 8.0
    if n < 500:   return 5.0
    if n < 2000:  return 3.5
    if n < 10000: return 2.0
    return 1.0


def load_sparse(path):
    A = mmread(str(path))
    return A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)


def with_diagonal(mat):
    mat = mat.tolil()
    for i in range(min(mat.shape)):
        mat[i, i] = 1
    return mat.tocsr()


def count_fbm(mat):
    coo  = mat.tocoo()
    mask = coo.col > coo.row
    return int(mask.sum()), int((coo.col[mask] - coo.row[mask]).sum())


def resolve_mtx(stem, repo_root):
    data_dir = repo_root / "data"
    exact = data_dir / (stem + ".mtx")
    if exact.exists():
        return exact
    for f in data_dir.glob("*.mtx"):
        if f.stem.lower() == stem.lower():
            return f
    raise SystemExit(f"Matrix '{stem}' not found in {data_dir}")


def find_permuted(matrix_path, repo_root, perm_arg):
    """Resolve the permuted matrix path for a given matrix."""
    if perm_arg:
        return Path(perm_arg)
    candidate = repo_root / "out" / "permuted" / (matrix_path.stem + "_permuted.mtx")
    if candidate.exists():
        return candidate
    return None


# ── drawing ───────────────────────────────────────────────────────────────────

def draw_scc_boxes(ax, block_sizes,
                   facecolor="tab:orange", edgecolor="black",
                   alpha=0.25, lw=0.8):
    pos = 0
    for size in block_sizes:
        rect = patches.Rectangle(
            (pos - 0.5, pos - 0.5), size, size,
            linewidth=lw,
            edgecolor=edgecolor,
            facecolor=facecolor,
            alpha=alpha,
            zorder=2
        )
        ax.add_patch(rect)
        pos += size


def add_stats_annotation(ax, fbm, tfbd):
    text = f"|FBM| = {fbm}\nTFBD = {tfbd}"
    ax.text(
        0.98, 0.98, text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="right",
        multialignment="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="lightgrey", alpha=0.9)
    )


def plot_spy(ax, mat, markersize=2, block_sizes=None,
             fbm=None, tfbd=None, show_ylabel=True):
    mat = with_diagonal(mat)
    if block_sizes is not None:
        draw_scc_boxes(ax, block_sizes)
    ax.spy(mat, markersize=markersize, color="steelblue", zorder=3)
    ax.set_xlabel("columns", fontsize=9, fontweight="bold")
    if show_ylabel:
        ax.set_ylabel("rows", fontsize=9, fontweight="bold")
    else:
        ax.set_ylabel("")
        ax.tick_params(axis="y", left=True, labelleft=False)
    ax.tick_params(labelsize=7)
    if fbm is not None and tfbd is not None:
        add_stats_annotation(ax, fbm, tfbd)


def add_matrix_title(fig, ax, title):
    pos = ax.get_position()
    x = pos.x0 - 0.075
    y = pos.y0 + pos.height / 2
    fig.text(x, y, title,
             fontsize=12, fontweight="bold",
             ha="center", va="center",
             rotation=90)


def load_permuted_matrix(perm_path, A, n, perm_fmt, one_based):
    """Load the permuted matrix from either a .mtx file or a permutation vector."""
    perm_path = Path(perm_path)
    if not perm_path.exists():
        return None
    try:
        if perm_path.suffix.lower() == ".mtx":
            B = load_sparse(perm_path)
            if B.shape != (n, n):
                raise SystemExit(f"Permuted matrix shape {B.shape} doesn't match original ({n},{n})")
            return B
        else:
            perm = load_permutation(str(perm_path), n, fmt=perm_fmt, one_based=one_based)
            return A[perm, :][:, perm]
    except Exception:
        perm = load_permutation(str(perm_path), n, fmt=perm_fmt, one_based=one_based)
        return A[perm, :][:, perm]


# ── single matrix mode ────────────────────────────────────────────────────────

def run_single(args, repo_root):
    matrix_path = Path(args.matrix)
    if not matrix_path.exists():
        alt = repo_root / "data" / args.matrix
        if alt.exists():
            matrix_path = alt
        else:
            raise SystemExit(f"Matrix file not found: {args.matrix}")

    stem      = matrix_path.stem
    sccs_path = repo_root / "out" / "sccs" / (stem + "_sccs.txt")

    A = load_sparse(matrix_path)
    n_rows, n_cols = A.shape
    if n_rows != n_cols:
        raise SystemExit("This script expects a square matrix for symmetric reordering (P^T A P).")
    n  = n_rows
    ms = args.markersize if args.markersize != 1.0 else auto_markersize(n)

    block_sizes = load_scc_blocks(sccs_path) if sccs_path.exists() else None
    fbm0, tfbd0 = count_fbm(A)

    perm_path = find_permuted(matrix_path, repo_root, args.perm)
    A_perm    = load_permuted_matrix(perm_path, A, n, args.perm_format, args.one_based) if perm_path else None

    if A_perm is not None:
        fbm1, tfbd1 = count_fbm(A_perm)
        fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2), squeeze=False)
        axes = axes[0]
        axes[0].set_title("Original",    fontsize=15, fontweight="bold", pad=13)
        axes[1].set_title("Partitioned", fontsize=15, fontweight="bold", pad=13)
        plt.tight_layout(w_pad=2.0)
        plot_spy(axes[0], A, markersize=ms, fbm=fbm0, tfbd=tfbd0, show_ylabel=True)
        plot_spy(axes[1], A_perm, markersize=ms, block_sizes=block_sizes,
                 fbm=fbm1, tfbd=tfbd1, show_ylabel=False)
        add_matrix_title(fig, axes[0], stem)
    else:
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.set_title("Original", fontsize=15, fontweight="bold", pad=13)
        plt.tight_layout()
        plot_spy(ax, A, markersize=ms, fbm=fbm0, tfbd=tfbd0, show_ylabel=True)
        add_matrix_title(fig, ax, stem)

    if args.output:
        fig.savefig(args.output, dpi=200, bbox_inches="tight")
        print(f"Saved figure to {args.output}")
    else:
        plt.show()


# ── multi-matrix comparison mode ─────────────────────────────────────────────

def run_compare(args, repo_root):
    output_dir = Path(args.output_dir) if args.output_dir else repo_root / "out" / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, stem in enumerate(args.compare):
        ms               = args.markersizes[i] if (args.markersizes and i < len(args.markersizes)) else None
        show_col_headers = (i == 0)

        mtx_path  = resolve_mtx(stem, repo_root)
        perm_path = repo_root / "out" / "permuted" / (mtx_path.stem + "_permuted.mtx")
        sccs_path = repo_root / "out" / "sccs"     / (mtx_path.stem + "_sccs.txt")

        A  = load_sparse(mtx_path)
        n  = A.shape[0]
        ms = ms if ms is not None else auto_markersize(n)

        block_sizes = load_scc_blocks(sccs_path) if sccs_path.exists() else None
        fbm0, tfbd0 = count_fbm(A)

        fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), squeeze=False)
        axes = axes[0]

        if show_col_headers:
            axes[0].set_title("Original",    fontsize=15, fontweight="bold", pad=13)
            axes[1].set_title("Partitioned", fontsize=15, fontweight="bold", pad=13)

        plt.tight_layout(h_pad=3.0, w_pad=2.0)

        plot_spy(axes[0], A, markersize=ms, fbm=fbm0, tfbd=tfbd0, show_ylabel=True)

        if perm_path.exists():
            A_perm      = load_sparse(perm_path)
            fbm1, tfbd1 = count_fbm(A_perm)
            plot_spy(axes[1], A_perm, markersize=ms, block_sizes=block_sizes,
                     fbm=fbm1, tfbd=tfbd1, show_ylabel=False)
        else:
            axes[1].set_visible(False)
            print(f"  Warning: no permuted matrix found for {stem}")

        add_matrix_title(fig, axes[0], mtx_path.stem)

        out_path = output_dir / f"{mtx_path.stem}.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved → {out_path}")


# ── entry point ───────────────────────────────────────────────────────────────

def main():
    repo_root = Path(__file__).resolve().parent.parent

    p = argparse.ArgumentParser(description="Plot sparse matrix before/after permutation")

    # single matrix mode
    p.add_argument("matrix", nargs="?", default=None,
                   help="Path to .mtx matrix file (single matrix mode)")
    p.add_argument("--perm",
                   help="Path to permutation file (npy, plain text, or .mtx)")
    p.add_argument("--perm-format", choices=["new2old", "old2new"], default="new2old",
                   help="Format of permutation file (default: new2old)")
    p.add_argument("--one-based", action="store_true",
                   help="Treat permutation as 1-based indices")
    p.add_argument("--output",
                   help="Output image file (png, pdf). If omitted, show interactively")
    p.add_argument("--markersize", type=float, default=1.0,
                   help="Marker size for spy plots (default: auto-scaled)")

    # multi-matrix comparison mode
    p.add_argument("--compare", nargs="+", metavar="STEM",
                   help="Matrix stems for multi-matrix comparison mode")
    p.add_argument("--markersizes", nargs="+", type=float, metavar="SIZE",
                   help="Per-matrix marker sizes for --compare mode")
    p.add_argument("--output-dir", default=None,
                   help="Output directory for PNGs in --compare mode (default: out/plots/)")

    args = p.parse_args()

    if args.compare:
        run_compare(args, repo_root)
    elif args.matrix:
        run_single(args, repo_root)
    else:
        p.print_help()


if __name__ == "__main__":
    main()