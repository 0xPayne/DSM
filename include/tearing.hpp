#pragma once
#include "sparse.hpp"
#include <vector>

namespace SparseLib::Tearing {

    std::vector<std::vector<int>> tearAllSCCs(
        const std::vector<std::vector<int>>& sccs,
        const Sparse::CSCMatrix& matrix
    );

}
