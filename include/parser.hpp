#pragma once
#include "sparse.hpp"
#include <string>

namespace SparseLib {
namespace IO {
    // Loads binary DSM data from a matrix market file 
    // and constructs a valueless CSC matrix.
    Sparse::CSCMatrix loadFromFile(const std::string& filepath);
    // Save a valueless CSC matrix to a Matrix Market (.mtx) file
    void saveToFile(const Sparse::CSCMatrix& matrix, const std::string& filepath);
}}
