// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include "sparse.hpp"
#include <vector>

namespace SparseLib::SCC {

    // find SCCs using Tarjan's or Kosaraju's algorithm
    std::vector<std::vector<int>> tarjanSCC(const Sparse::CSCMatrix& matrix);
    std::vector<std::vector<int>> kosarajuSCC(const Sparse::CSCMatrix& matrix);

    // build the DAG of SCCs
    Sparse::CSCMatrix condensationGraph(const Sparse::CSCMatrix& matrix, const std::vector<std::vector<int>>& sccs);

    // topological sort (Kahn's algorithm)
    std::vector<int> topologicalSort(const Sparse::CSCMatrix& matrix);

} // namespace SparseLib::SCC
