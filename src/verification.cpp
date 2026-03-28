// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/verification.hpp"
#include <algorithm>
#include <utility>

namespace SparseLib::Verify {

Result verifyPermutation(
    const Sparse::CSCMatrix& original,
    const Sparse::CSCMatrix& permuted,
    const std::vector<int>& perm)
{
    const int n = original.cols;

    // 1. perm must be a valid permutation of [0, n)
    if (static_cast<int>(perm.size()) != n)
        return {false, "perm size (" + std::to_string(perm.size()) +
                        ") != matrix dimension (" + std::to_string(n) + ")"};

    std::vector<bool> seen(n, false);
    for (int v : perm) {
        if (v < 0 || v >= n)
            return {false, "perm contains out-of-range value " + std::to_string(v)};
        if (seen[v])
            return {false, "perm contains duplicate value " + std::to_string(v)};
        seen[v] = true;
    }

    // 2. Dimensions and nnz must match
    if (permuted.rows != original.rows || permuted.cols != original.cols)
        return {false, "dimension mismatch: original " +
                        std::to_string(original.rows) + "x" + std::to_string(original.cols) +
                        ", permuted " +
                        std::to_string(permuted.rows) + "x" + std::to_string(permuted.cols)};

    if (permuted.nnz != original.nnz)
        return {false, "nnz mismatch: original " + std::to_string(original.nnz) +
                        ", permuted " + std::to_string(permuted.nnz)};

    // 3. edge sets must match under the permutation
    using Edge = std::pair<int,int>;

    auto collectEdges = [](const Sparse::CSCMatrix& m) {
        std::vector<Edge> edges;
        edges.reserve(m.nnz);
        for (int col = 0; col < m.cols; col++)
            for (int p = m.col_ptrs[col]; p < m.col_ptrs[col + 1]; p++)
                edges.emplace_back(m.row_indices[p], col);
        return edges;
    };

    std::vector<Edge> orig_edges = collectEdges(original);

    // un-permute back to original space
    std::vector<Edge> perm_edges;
    perm_edges.reserve(permuted.nnz);
    for (int col = 0; col < permuted.cols; col++)
        for (int p = permuted.col_ptrs[col]; p < permuted.col_ptrs[col + 1]; p++)
            perm_edges.emplace_back(perm[permuted.row_indices[p]], perm[col]);

    std::sort(orig_edges.begin(), orig_edges.end());
    std::sort(perm_edges.begin(), perm_edges.end());

    if (orig_edges != perm_edges) {
        size_t i = 0;
        while (i < orig_edges.size() && i < perm_edges.size() && orig_edges[i] == perm_edges[i])
            i++;
        std::string detail;
        if (i < orig_edges.size() && i < perm_edges.size())
            detail = "first mismatch at index " + std::to_string(i) +
                     ": original (" + std::to_string(orig_edges[i].first) + "," +
                     std::to_string(orig_edges[i].second) + ") vs permuted (" +
                     std::to_string(perm_edges[i].first) + "," +
                     std::to_string(perm_edges[i].second) + ")";
        else
            detail = "edge count diverged at index " + std::to_string(i);
        return {false, "edge sets differ: " + detail};
    }

    return {true, "ok"};
}

} // namespace SparseLib::Verify
