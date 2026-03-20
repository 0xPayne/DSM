#pragma once
#include "sparse.hpp"

#include <string>
#include <vector>

namespace SparseLib::Verify {

    struct Result {
        bool passed;
        std::string message;
    };

    // Verify that `permuted` is exactly the symmetric permutation P^T A 
    // of `original`, with no edges lost or introduced.
    Result verifyPermutation(
        const Sparse::CSCMatrix& original,
        const Sparse::CSCMatrix& permuted,
        const std::vector<int>& perm);

}
