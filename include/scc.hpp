#pragma once
#include "sparse.hpp"
#include <vector>

namespace SparseLib {

namespace SCC {

    // Returns a list of strongly connected components; each component is a vector of vertex indices
    // CSC matrix is interpreted as outgoing edges of vertex u are in column u
    std::vector<std::vector<int>> tarjanSCC(const Sparse::CSCMatrix& matrix);

    // Build the condensation graph of the Strongly Connected Components
    Sparse::CSCMatrix condensationGraph(const Sparse::CSCMatrix& matrix, const std::vector<std::vector<int>>& sccs);

}} 
