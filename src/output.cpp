// CPSC 482 : DSM Project 
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/output.hpp"
#include "../include/parser.hpp"
#include "../include/ux.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

namespace SparseLib::Output {

void savePermutedMatrix(const std::string &execPath,
                        const std::string &inputPath,
                        const Sparse::CSCMatrix &permuted) {
    try {
        fs::path outDir = SparseLib::UX::projectRoot(execPath) / "out" / "permuted";
        fs::create_directories(outDir);
        fs::path outPath =
        outDir / (fs::path(inputPath).stem().string() + "_permuted.mtx");
        SparseLib::IO::saveToFile(permuted, outPath.string());
        std::cout << "Wrote permuted matrix to: " << outPath.string() << "\n";
    } catch (const std::exception &e) {
    std::cerr << "Warning: failed to write permuted matrix: " << e.what()
              << "\n";
    }
}

void saveSCCBlocks(const std::string &execPath,
                   const std::string &inputPath,
                   const std::vector<std::vector<int>> &sccs,
                   const std::vector<int> &topoOrder) {
    try {
        fs::path outDir = SparseLib::UX::projectRoot(execPath) / "out" / "sccs";
        fs::create_directories(outDir);
        fs::path outPath = outDir / (fs::path(inputPath).stem().string() + "_sccs.txt");

        std::ofstream f(outPath);
        if (!f)
            throw std::runtime_error("Cannot open " + outPath.string());

        // One SCC size per line, in topological order
        for (int idx : topoOrder)
            f << sccs[idx].size() << "\n";

        std::cout << "Wrote SCC blocks to: " << outPath.string() << "\n";
    } catch (const std::exception &e) {
        std::cerr << "Warning: failed to write SCC blocks: " << e.what() << "\n";
    }
}

}