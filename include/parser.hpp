// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include "sparse.hpp"
#include <string>

namespace SparseLib::IO {

    // Loads binary DSM data from a matrix market file
    // and constructs a valueless CSC matrix.
    Sparse::CSCMatrix loadFromFile(const std::string& filepath);

    // Save a valueless CSC matrix to a Matrix Market (.mtx) file
    void saveToFile(const Sparse::CSCMatrix& matrix, const std::string& filepath);

} // namespace SparseLib::IO
