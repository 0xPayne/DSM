#include "../include/parser.hpp"


namespace SparseLib::IO {
    Sparse::CSCMatrix loadFromFile(const std::string& filepath) {
        int rows, cols, nnz;

        // Using pair<int, int> where .first = row, .second = col
        std::vector<std::pair<int, int>> duplets;


        std::ifstream file(filepath);
        if (!file.is_open()) throw std::runtime_error("Could not open file: " + filepath);

        std::string line;
        std::getline(file, line);

        // Skip comments
        while (std::getline(file, line) && line[0] == '%');

        // Read dimensions
        std::stringstream ss(line);
        if (!(ss >> rows >> cols >> nnz)) {
            throw std::runtime_error("Invalid Matrix Market header in: " + filepath);
        }

        duplets.resize(nnz);
        for (int i = 0; i < nnz; i++) {
            int r, c;
            // Matrix Market files often have a 3rd column for values; 
            // since this is valueless, we just read r and c and ignore the rest of the line.
            if (!(file >> r >> c)) throw std::runtime_error("Read error at line " + std::to_string(i));
            duplets[i] = {r - 1, c - 1}; // Convert to 0-based indexing
            
            // Consume the rest of the line (in case there are weights/values)
            std::string dummy;
            std::getline(file, dummy);
        }

        std::sort(duplets.begin(), duplets.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            if (a.second != b.second) return a.second < b.second;
            return a.first < b.first;
        });

        Sparse::CSCMatrix matrix{rows, cols, nnz};
        matrix.row_indices.resize(nnz);
        matrix.col_ptrs.assign(cols + 1, 0);


        // Fill row indices and build the histogram for column pointers
        for (int i = 0; i < nnz; i++) {
            matrix.row_indices[i] = duplets[i].first;  // Row
            matrix.col_ptrs[duplets[i].second + 1]++;   // Column histogram
        }

        // Prefix sum
        for (int i = 0; i < cols; i++) {
            matrix.col_ptrs[i + 1] += matrix.col_ptrs[i];
        }

        return matrix;
    }
}

