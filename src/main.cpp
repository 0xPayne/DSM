#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

struct CSRMatrix {
    int rows, cols, nnz;
    std::vector<int> row_ptr;    // Size: rows + 1
    std::vector<int> col_ind;    // Size: nnz
    std::vector<double> values;  // Size: nnz (if DSM has weights)
};

struct Entry { 
    int r, c; 
    double v; 
};

CSRMatrix loadMatrixMarket(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string line;
    // 1. Skip comments and header
    while (std::getline(file, line) && line[0] == '%');

    // 2. Read dimensions
    int num_rows, num_cols, num_nnz;
    std::stringstream ss(line);
    ss >> num_rows >> num_cols >> num_nnz;

    CSRMatrix matrix;
    matrix.rows = num_rows;
    matrix.cols = num_cols;
    matrix.nnz = num_nnz;
    
    // We'll use a temporary count array for the rows
    std::vector<int> counts(num_rows, 0);
    
    // Store coordinates temporarily to avoid a third pass
    // Note: If memory is EXTREMELY tight, do a second pass of the FILE instead
    std::vector<Entry> entries(num_nnz);

    for (int i = 0; i < num_nnz; ++i) {
        file >> entries[i].r >> entries[i].c >> entries[i].v;
        entries[i].r--; // Convert 1-based to 0-based
        entries[i].c--;
        counts[entries[i].r]++;
    }

    // 3. Build row_ptr
    matrix.row_ptr.resize(num_rows + 1, 0);
    for (int i = 0; i < num_rows; ++i) {
        matrix.row_ptr[i + 1] = matrix.row_ptr[i] + counts[i];
    }

    // 4. Fill col_ind and values
    matrix.col_ind.resize(num_nnz);
    matrix.values.resize(num_nnz);
    
    // Use a copy of row_ptr to track where to insert next element in each row
    std::vector<int> current_pos = matrix.row_ptr;

    for (int i = 0; i < num_nnz; ++i) {
        int row = entries[i].r;
        int dest = current_pos[row];
        matrix.col_ind[dest] = entries[i].c;
        matrix.values[dest] = entries[i].v;
        current_pos[row]++;
    }

    return matrix;
}
void compareMemoryUsage(long long rows, long long cols, long long nnz) {
    // Normal 2D Array: Rows * Cols * sizeof(double)
    // We use double (8 bytes) as per your Entry struct
    double normalBytes = (double)rows * cols * sizeof(double);
    double normalGB = normalBytes / (1024 * 1024 * 1024);

    // CSR Storage:
    // row_ptr: (rows + 1) * sizeof(int)
    // col_ind: nnz * sizeof(int)
    // values:  nnz * sizeof(double)
    double csrBytes = ((rows + 1) * sizeof(int)) + (nnz * sizeof(int)) + (nnz * sizeof(double));
    double csrMB = csrBytes / (1024 * 1024);

    std::cout << "\n=== MEMORY CONSUMPTION COMPARISON ===" << std::endl;
    std::cout << "Matrix Size: " << rows << " x " << cols << std::endl;
    std::cout << "Non-zero Entries: " << nnz << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "NORMAL (2D Array) Storage:  ~" << normalGB << " GB" << std::endl;
    std::cout << "CSR Storage:                ~" << csrMB << " MB" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "CSR is approximately " << (normalGB * 1024) / csrMB << "x more efficient." << std::endl;
    std::cout << "=====================================\n" << std::endl;
}
void printMetrics(const CSRMatrix& mat) {
    long long feedback = 0;
    long long dist_sum = 0;
    for (int i = 0; i < mat.rows; i++) {
        for (int j = mat.row_ptr[i]; j < mat.row_ptr[i+1]; j++) {
            if (mat.col_ind[j] > i) {
                feedback++;
                dist_sum += (mat.col_ind[j] - i);
            }
        }
    }
    std::cout << "--- Initial Metrics ---" << std::endl;
    std::cout << "Non-zeros above diagonal: " << feedback << std::endl;
    std::cout << "Sum of distances: " << dist_sum << std::endl;
}

int main() {
    try {
        std::string filename = "data/web-BerkStan/web-BerkStan.mtx";
        
        // We read only the header first to show memory comparison before loading
        std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line) && line[0] == '%');
        int r, c, n;
        std::stringstream ss(line);
        ss >> r >> c >> n;
        file.close();

        // 1. Show the Comparison
        compareMemoryUsage(r, c, n);

        // 2. Load with CSR
        auto start = std::chrono::high_resolution_clock::now();
        CSRMatrix matrix = loadMatrixMarket(filename);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Loaded into CSR in: " << std::chrono::duration<double>(end-start).count() << "s\n";

        // 3. Metrics and SCC
        printMetrics(matrix);
        // ... (remaining SCC code)

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}