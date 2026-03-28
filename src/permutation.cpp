// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/permutation.hpp"
#include <algorithm>

namespace SparseLib::Permutation {

std::vector<int> buildPermutation(
    const std::vector<std::vector<int>>& sccs,
    const std::vector<int>& topo_order)
{
    std::vector<int> perm;

    // append SCC vertices in topological order
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

    // remap entries via histogram, prefix sum, scatter

    Sparse::CSCMatrix result(n, n, matrix.nnz);
    result.col_ptrs.assign(n + 1, 0);
    result.row_indices.resize(matrix.nnz);

    // histogram: count entries per new column
    for (int old_col = 0; old_col < n; ++old_col) {
        int count   = matrix.col_ptrs[old_col + 1] - matrix.col_ptrs[old_col];
        int new_col = inv[old_col];
        result.col_ptrs[new_col + 1] += count;
    }

    // prefix sum
    for (int i = 0; i < n; ++i) {
        result.col_ptrs[i + 1] += result.col_ptrs[i];
    }

    // scatter entries into new positions
    std::vector<int> write_pos(result.col_ptrs.begin(), result.col_ptrs.begin() + n);

    for (int old_col = 0; old_col < n; ++old_col) {
        int new_col = inv[old_col];
        for (int p = matrix.col_ptrs[old_col]; p < matrix.col_ptrs[old_col + 1]; ++p) {
            int new_row = inv[matrix.row_indices[p]];
            result.row_indices[write_pos[new_col]++] = new_row;
        }
    }

    // sort rows within each column
    for (int col = 0; col < n; ++col) {
        std::sort(result.row_indices.begin() + result.col_ptrs[col],
                  result.row_indices.begin() + result.col_ptrs[col + 1]);
    }

    return result;
}

} // namespace SparseLib::Permutation
