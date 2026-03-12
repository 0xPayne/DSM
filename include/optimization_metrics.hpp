#pragma once
#include "sparse.hpp"

namespace SparseLib {
namespace Metrics {

    // Calculates the Feedback Marks (FBM) count.
    int countFBM(const Sparse::CSCMatrix& dsm);

    // Distance of FBM from the diagonal sum.
    double fbmDiagonalDistance(const Sparse::CSCMatrix& dsm);
}}
