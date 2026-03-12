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

    struct CSRMatrix {
        int rows;
        int cols;
        int nnz;
        std::vector<int> row_ptrs;
        std::vector<int> col_indices;
        CSRMatrix(int r = 0, int c = 0, int n = 0) : rows(r), cols(c), nnz(n) {}
    };

    CSRMatrix transpose(CSCMatrix&& mat);
    CSCMatrix transpose(CSRMatrix&& mat);

    CSCMatrix add(const CSCMatrix& a, const CSCMatrix& b);
    CSRMatrix add(const CSRMatrix& a, const CSRMatrix& b);

    CSCMatrix multiply(const CSCMatrix& a, const CSCMatrix& b);
    CSRMatrix multiply(const CSRMatrix& a, const CSRMatrix& b);

}
