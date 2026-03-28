// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include <vector>

namespace Sparse {

    // Sparse Storage Structures : Compressed Sparse Column & Compressed Sparse Row.
    struct CSCMatrix {
        int rows;
        int cols;
        int nnz;
        std::vector<int> col_ptrs;
        std::vector<int> row_indices;
        CSCMatrix(int r = 0, int c = 0, int n = 0) : rows(r), cols(c), nnz(n) {}
    };

    CSCMatrix transpose(CSCMatrix&& mat);

} // namespace Sparse
