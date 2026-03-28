// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/sparse.hpp"

namespace Sparse {

Sparse::CSCMatrix transpose(Sparse::CSCMatrix&& mat) {
  const int n = mat.cols;
  Sparse::CSCMatrix result(n, n, mat.nnz);
  result.col_ptrs.assign(n + 1, 0);
  result.row_indices.assign(mat.nnz, 0);

  for (int i = 0; i < mat.nnz; ++i) {
    result.col_ptrs[mat.row_indices[i] + 1]++;
  }

  for (int i = 0; i < n; ++i) {
    result.col_ptrs[i + 1] += result.col_ptrs[i];
  }

  std::vector<int> current_ptrs = result.col_ptrs;
  for (int u = 0; u < n; ++u) {
    for (int p = mat.col_ptrs[u]; p < mat.col_ptrs[u + 1]; ++p) {
      int v = mat.row_indices[p];
      result.row_indices[current_ptrs[v]++] = u;
    }
  }

  return result;
}

} // namespace Sparse
