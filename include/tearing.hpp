// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include "sparse.hpp"
#include <vector>

namespace SparseLib::Tearing {

    std::vector<std::vector<int>> tearAllSCCs(
        const std::vector<std::vector<int>>& sccs,
        const Sparse::CSCMatrix& matrix
    );

}
