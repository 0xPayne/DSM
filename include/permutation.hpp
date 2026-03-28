// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include "sparse.hpp"
#include <vector>

namespace SparseLib::Permutation {

    // build permutation from SCC decomposition + topo order
    std::vector<int> buildPermutation(
        const std::vector<std::vector<int>>& sccs,
        const std::vector<int>& topo_order
    );

    // invert a permutation
    std::vector<int> invertPermutation(const std::vector<int>& perm);

    // apply symmetric permutation P^T A P
    Sparse::CSCMatrix applyPermutation(
        const Sparse::CSCMatrix& matrix,
        const std::vector<int>& perm
    );

}
