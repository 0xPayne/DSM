#pragma once
#include "sparse.hpp"
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

namespace SparseLib {
namespace IO {
    // Loads binary DSM data from a matrix market file 
    // and constructs a valueless CSC matrix.
    Sparse::CSCMatrix loadFromFile(const std::string& filepath);
}}
