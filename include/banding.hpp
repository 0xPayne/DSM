// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include "sparse.hpp"
#include <vector>

namespace SparseLib::Banding {

    // Two-phase banding per SCC: RCM reordering + adjacent-swap hill climbing.
    std::vector<std::vector<int>> bandAllSCCs(
        const std::vector<std::vector<int>>& torn_sccs,
        const Sparse::CSCMatrix& matrix
    );

}
