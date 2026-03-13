#include "../include/optimization_metrics.hpp"

namespace SparseLib::Metrics {

    int countFBM(const Sparse::CSCMatrix& dsm) {
	int c = 0;
	for (int col = 0; col < dsm.cols; col++) {
	    int start = dsm.col_ptrs[col];
	    int end = dsm.col_ptrs[col + 1];
	    for (int idx = start; idx < end; idx++) {
		int row = dsm.row_indices[idx];
		if (row < col) {
		    c++;
		}
	    }
	}
	return c;
    }

    long long fbmDiagonalDistance(const Sparse::CSCMatrix& dsm) {
    long long total = 0;
    for (int col = 0; col < dsm.cols; col++) {
        const int start = dsm.col_ptrs[col];
        const int end   = dsm.col_ptrs[col + 1];
        for (int idx = start; idx < end; idx++) {
            const int row = dsm.row_indices[idx];
            if (row < col) {
                total += col - row;
            }
        }
    }
    return total;
}
}

