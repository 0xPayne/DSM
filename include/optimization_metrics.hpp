// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include "sparse.hpp"

namespace SparseLib::Metrics {

    // count of above-diagonal entries (feedback marks)
    int countFBM(const Sparse::CSCMatrix& dsm);

    // sum of distances of feedback marks from diagonal
    long long fbmDiagonalDistance(const Sparse::CSCMatrix& dsm);

} // namespace SparseLib::Metrics
