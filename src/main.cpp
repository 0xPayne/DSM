#include <iostream>
#include "../include/parser.hpp"
#include "../include/optimization_metrics.hpp"

int main() {
    
    // To Do:
    // 1. Read in MM file with the parser.
    // 2. Analyze its stock metrics.
    // 3. Extract SCC's ..........
    std::string filepath = "data/easy-example.mtx";

    std::cout << "Processing (CSR) " << filepath << "..." << std::endl;

    Sparse::CSCMatrix cscMatrix = SparseLib::IO::loadFromFile(filepath);

    std::cout << "# FBM = " << SparseLib::Metrics::countFBM(cscMatrix) << std::endl;
    std::cout << "Sum Distance of DBM = " << SparseLib::Metrics::fbmDiagonalDistance(cscMatrix) << std::endl;
    

    return 0;
}
