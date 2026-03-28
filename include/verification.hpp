// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include "sparse.hpp"

#include <string>
#include <vector>

namespace SparseLib::Verify {

    struct Result {
        bool passed;
        std::string message;
    };

    // check that permuted == P^T A P with no lost or extra edges
    Result verifyPermutation(
        const Sparse::CSCMatrix& original,
        const Sparse::CSCMatrix& permuted,
        const std::vector<int>& perm);

}
