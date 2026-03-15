#pragma once
#include "sparse.hpp"
#include <vector>

namespace SparseLib::Permutation {

    // Build a permutation vector from SCC decomposition + topological ordering.
    //
    // sccs[i]       = list of vertices in SCC i
    // topo_order[k] = which SCC index should appear at position k
    //
    // Returns perm where perm[new_position] = old_vertex.
    std::vector<int> buildPermutation(
        const std::vector<std::vector<int>>& sccs,
        const std::vector<int>& topo_order
    );

    // Invert a permutation: if perm[new] = old, returns inv where inv[old] = new.
    std::vector<int> invertPermutation(const std::vector<int>& perm);

    // Apply symmetric permutation P^T A P to a CSC matrix.
    // Reorders both rows and columns so the matrix reflects the new vertex ordering.
    Sparse::CSCMatrix applyPermutation(
        const Sparse::CSCMatrix& matrix,
        const std::vector<int>& perm
    );

}
