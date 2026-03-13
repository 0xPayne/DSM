#include <iostream>
#include "../include/parser.hpp"
#include "../include/scc.hpp"
#include "../include/optimization_metrics.hpp"

int main() {
    
    // To Do:
    // 1. Read in MM file with the parser.
    // 2. Analyze its stock metrics.
    // 3. Extract SCC's ..........
    std::string filepath = "data/Tina_AskCal.mtx";

    std::cout << "Processing (CSR) " << filepath << "..." << std::endl;

    Sparse::CSCMatrix cscMatrix = SparseLib::IO::loadFromFile(filepath);

    std::vector<std::vector<int>> sccs = SparseLib::SCC::tarjanSCC(cscMatrix);

    std::cout << "Number of SCCs: " << sccs.size() << std::endl;

    return 0;
}
