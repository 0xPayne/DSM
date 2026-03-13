#include <iostream>
#include "../include/parser.hpp"
#include "../include/scc.hpp"
#include "../include/optimization_metrics.hpp"

int main() {
    
    // To Do:
    // 1. Read in MM file with the parser.
    // 2. Analyze its stock metrics.
    // 3. Extract SCC's 
    // 4. Create condensation graph
    // 5. Perform topologisation graph cal sort on conden

    std::string filepath = "data/HEP-th-new.mtx";

    std::cout << "Processing (CSR) " << filepath << "..." << std::endl;

    Sparse::CSCMatrix cscMatrix = SparseLib::IO::loadFromFile(filepath);

    std::vector<std::vector<int>> sccs = SparseLib::SCC::tarjanSCC(cscMatrix);

    int numSccs = sccs.size();

    std::cout << "Number of SCCs: " << numSccs << std::endl;

    Sparse::CSCMatrix condense = SparseLib::SCC::condensationGraph(cscMatrix, sccs);

    return 0;
}
