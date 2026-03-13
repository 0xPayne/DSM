#include "../include/parser.hpp"
#include "../include/scc.hpp"
#include "../include/ux.hpp"
#include "../include/optimization_metrics.hpp"

#include <iostream>

int main(int argc, char* argv[]) {
    try {
        // 1. Select a matrix market file and parse it into compressed row/column.
        std::string filepath = SparseLib::UX::selectFileFromDirectory(argv[0]);

        std::cout << "\nProcessing (CSC) " << filepath << "..." << std::endl;
        Sparse::CSCMatrix dsm = SparseLib::IO::loadFromFile(filepath);
        
        // 2. Analyze stock metrics: Number of Feedback marks, Sum of Distance of Feedback marks from the diagonal.
        std::cout << "Number of Feedback Marks: " << SparseLib::Metrics::countFBM(dsm) << std::endl;
        std::cout << "Sum of Distances of FBM's from Diagonal: " << SparseLib::Metrics::fbmDiagonalDistance(dsm) << std::endl;

        // 3. Deduce the SCC's of the DSM.
        std::vector<std::vector<int>> sccs = SparseLib::SCC::tarjanSCC(dsm);
        std::cout << "Number of SCCs: " << sccs.size() << std::endl;

        // 4. Create condensation graph
        Sparse::CSCMatrix condense = SparseLib::SCC::condensationGraph(dsm, sccs);

        // 5. Perform topological sort on condensation graph 

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
