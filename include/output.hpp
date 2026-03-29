// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include <sparse.hpp>
#include <string>
#include <vector>

namespace SparseLib::Output {

// save permuted matrix to out/permuted/<matrix-name>_permuted.mtx
void savePermutedMatrix(const std::string &execPath,
                const std::string &inputPath,
                const Sparse::CSCMatrix &permuted
);

// save scc blocks to out/sccs/<matrix-name>_sccs.txt
// one integer per line, each being the size of one scc block
// in topological order, matching the rows and cols of the permuted matrix
void saveSCCBlocks(const std::string &execPath,
                const std::string &inputPath,
                const std::vector<std::vector<int>> &sccs,
                const std::vector<int> &topoOrder
);

} // namespace SparseLib::Output
