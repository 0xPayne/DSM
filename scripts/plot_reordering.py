#!/usr/bin/env python3
"""
Plot a sparse Matrix Market (.mtx) matrix before and after applying a permutation.

Usage examples:
  python scripts/plot_reordering.py data/easy-example.mtx
  python scripts/plot_reordering.py data/easy-example.mtx --perm perm.npy --output fig.png
  python scripts/plot_reordering.py data/easy-example.mtx --perm perm.txt --perm-format old2new

Permutation formats supported:
- "new2old" (default): array p of length n where p[new_index] = old_index (0-based)
- "old2new": array q where q[old_index] = new_index (0-based); the script will convert with np.argsort

If your permutation file is 1-based (Matlab-style), pass `--one-based`.
"""

import argparse
import numpy as np
import scipy.sparse as sp
from scipy.io import mmread
import matplotlib.pyplot as plt
import os
from pathlib import Path


def load_permutation(path, n, fmt="new2old", one_based=False):
    if path.endswith(".npy"):
        perm = np.load(path)
    else:
        perm = np.loadtxt(path, dtype=int)

    perm = perm.astype(int)
    if one_based:
        perm = perm - 1

    if perm.size != n:
        raise ValueError(f"Permutation length {perm.size} != matrix size {n}")

    if fmt == "old2new":
        # perm[old] = new  -> convert to new2old: new2old[new] = old
        perm = np.argsort(perm)
    elif fmt != "new2old":
        raise ValueError("Unknown perm format: choose 'new2old' or 'old2new'")

    return perm


def plot_matrix(ax, mat, title=None, markersize=1):
    # Use spy for large sparse matrices; convert to CSR for efficient row/col slicing
    ax.spy(mat, markersize=markersize)
    if title:
        ax.set_title(title)
    ax.set_xlabel('cols')
    ax.set_ylabel('rows')


def main():
    p = argparse.ArgumentParser(description="Plot sparse matrix before/after permutation")
    p.add_argument('matrix', help='Path to .mtx matrix file')
    p.add_argument('--perm', help='Path to permutation file (npy or plain text)')
    p.add_argument('--perm-format', choices=["new2old", "old2new"], default="new2old",
                   help="Format of provided permutation file")
    p.add_argument('--one-based', action='store_true', help='Treat permutation as 1-based indices')
    p.add_argument('--output', help='Output image file (png, pdf). If omitted, show interactively')
    p.add_argument('--markersize', type=float, default=1.0, help='Marker size for spy plots')
    args = p.parse_args()

    # Resolve matrix path: allow passing just the filename located in `data/`.
    repo_root = Path(__file__).resolve().parent.parent
    matrix_path = Path(args.matrix)
    if not matrix_path.exists():
        alt = repo_root / 'data' / args.matrix
        if alt.exists():
            matrix_path = alt
        else:
            raise SystemExit(f"Matrix file not found: {args.matrix}")
    # If no --perm provided, try to find a permuted matrix in out/permuted/
    if not args.perm:
        candidate = repo_root / 'out' / 'permuted' / (matrix_path.stem + '_permuted.mtx')
        if candidate.exists():
            args.perm = str(candidate)

    args.matrix = str(matrix_path)

    A = mmread(args.matrix)
    if not sp.issparse(A):
        A = sp.csr_matrix(A)
    else:
        A = A.tocsr()

    n_rows, n_cols = A.shape
    if n_rows != n_cols:
        raise SystemExit("This script expects a square matrix for symmetric reordering (P^T A P).")
    n = n_rows

    if args.perm:
        # If the provided --perm looks like a Matrix Market file, load it as a matrix.
        # This allows passing the C++-written `<name>_permuted.mtx` directly.
        try:
            if args.perm.lower().endswith('.mtx'):
                B = mmread(args.perm)
                if not sp.issparse(B):
                    B = sp.csr_matrix(B)
                B = B.tocsr()
                if B.shape != (n, n):
                    raise SystemExit(f"Permuted matrix shape {B.shape} doesn't match original {A.shape}")
                A_perm = B
            else:
                perm = load_permutation(args.perm, n, fmt=args.perm_format, one_based=args.one_based)
                # perm is new->old mapping: perm[new] = old
                A_perm = A[perm, :][:, perm]
        except Exception:
            # If mmread failed (file not a matrix), fall back to reading a permutation vector.
            perm = load_permutation(args.perm, n, fmt=args.perm_format, one_based=args.one_based)
            A_perm = A[perm, :][:, perm]
    else:
        perm = None
        A_perm = None

    # Plot side-by-side if perm provided, otherwise single plot
    if A_perm is not None:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        plot_matrix(axes[0], A, title='Original', markersize=args.markersize)
        plot_matrix(axes[1], A_perm, title='Permuted', markersize=args.markersize)
        plt.tight_layout()
    else:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        plot_matrix(ax, A, title='Original', markersize=args.markersize)
        plt.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=200)
        print(f"Saved figure to {args.output}")
    else:
        plt.show()


if __name__ == '__main__':
    main()
