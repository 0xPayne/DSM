#include "../include/permutation.hpp"
#include <algorithm>

namespace SparseLib::Permutation {

std::vector<int> buildPermutation(
    const std::vector<std::vector<int>>& sccs,
    const std::vector<int>& topo_order)
{
    std::vector<int> perm;

    // Walk SCCs in topological order; append each SCC's vertices.
    // This guarantees all inter-SCC edges become feed-forward (below diagonal).
    for (int scc_idx : topo_order) {
        for (int vertex : sccs[scc_idx]) {
            perm.push_back(vertex);
        }
    }

    return perm;
}

std::vector<int> invertPermutation(const std::vector<int>& perm) {
    int n = static_cast<int>(perm.size());
    std::vector<int> inv(n);

    // perm[new_pos] = old_vertex  -->  inv[old_vertex] = new_pos
    for (int new_pos = 0; new_pos < n; ++new_pos) {
        inv[perm[new_pos]] = new_pos;
    }

    return inv;
}

Sparse::CSCMatrix applyPermutation(
    const Sparse::CSCMatrix& matrix,
    const std::vector<int>& perm)
{
    const int n = matrix.cols;
    std::vector<int> inv = invertPermutation(perm);

    // We need to build a new CSC where every entry (old_row, old_col)
    // is remapped to (inv[old_row], inv[old_col]).
    //
    // Same 3-step pattern the parser uses:
    //   1. Histogram  -- count entries per new column
    //   2. Prefix sum -- convert counts to column pointers
    //   3. Scatter    -- place each entry at its correct position

    Sparse::CSCMatrix result(n, n, matrix.nnz);
    result.col_ptrs.assign(n + 1, 0);
    result.row_indices.resize(matrix.nnz);

    // Step 1: Each old column's entries all move to one new column.
    // So new_col = inv[old_col] gets exactly as many entries as old_col had.
    for (int old_col = 0; old_col < n; ++old_col) {
        int count   = matrix.col_ptrs[old_col + 1] - matrix.col_ptrs[old_col];
        int new_col = inv[old_col];
        result.col_ptrs[new_col + 1] += count;
    }

    // Step 2: Prefix sum to turn counts into pointers.
    for (int i = 0; i < n; ++i) {
        result.col_ptrs[i + 1] += result.col_ptrs[i];
    }

    // Step 3: Walk the original matrix and write each entry to its new location.
    // write_pos[c] tracks the next free slot in new column c.
    std::vector<int> write_pos(result.col_ptrs.begin(), result.col_ptrs.begin() + n);

    for (int old_col = 0; old_col < n; ++old_col) {
        int new_col = inv[old_col];
        for (int p = matrix.col_ptrs[old_col]; p < matrix.col_ptrs[old_col + 1]; ++p) {
            int new_row = inv[matrix.row_indices[p]];
            result.row_indices[write_pos[new_col]++] = new_row;
        }
    }

    // Sort row indices within each column for consistent ordering.
    for (int col = 0; col < n; ++col) {
        std::sort(result.row_indices.begin() + result.col_ptrs[col],
                  result.row_indices.begin() + result.col_ptrs[col + 1]);
    }

    return result;
}

} // namespace SparseLib::Permutation
