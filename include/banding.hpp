#pragma once
#include "sparse.hpp"
#include <vector>

namespace SparseLib::Banding {

    /// Apply bandwidth-reduction banding to each SCC to minimise TFBD
    /// (Total Feedback-mark Diagonal Distance) while preserving the FBM
    /// (Feedback Mark count) achieved by tearing.
    ///
    /// Two-phase refinement per non-trivial SCC:
    ///   Phase 1 -- Reverse Cuthill-McKee (RCM) reordering.
    ///   Phase 2 -- Adjacent-swap hill climbing on TFBD.
    ///
    /// torn_sccs : per-SCC vertex orderings produced by tearing.
    /// matrix    : the original (un-permuted) CSC adjacency matrix.
    /// Returns   : refined per-SCC vertex orderings.
    std::vector<std::vector<int>> bandAllSCCs(
        const std::vector<std::vector<int>>& torn_sccs,
        const Sparse::CSCMatrix& matrix
    );

}
