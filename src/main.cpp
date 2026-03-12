#include <iostream>
#include "parser.hpp"

int main() {
    
    // To Do:
    // 1. Read in MM file with the parser.
    // 2. Analyze its stock metrics.
    // 3. Extract SCC's ..........
    std::string filepath = "data/Tina_AskCal.mtx";

    std::cout << "Processing (CSR) " << filepath << "..." << std::endl;

    Sparse::CSCMatrix cscMatrix = SparseLib::IO::loadFromFile(filepath);
    

    return 0;
}
