#include "../include/benchmark.hpp"
#include "../include/optimization_metrics.hpp"
#include "../include/parser.hpp"
#include "../include/permutation.hpp"
#include "../include/scc.hpp"
#include "../include/ux.hpp"
#include "../include/verification.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Collect all .mtx files from the data/ directory
// ---------------------------------------------------------------------------
static std::vector<std::string> findMatrices(const std::string &execPath) {
  fs::path dataDir =
      fs::canonical(fs::path(execPath).parent_path() / "../data");
  std::vector<std::string> files;
  for (const auto &entry : fs::directory_iterator(dataDir))
    if (entry.path().extension() == ".mtx")
      files.push_back(entry.path().string());
  std::sort(files.begin(), files.end());
  return files;
}

// ---------------------------------------------------------------------------
// Interactive mode: original step-by-step flow, extended with permutation
// ---------------------------------------------------------------------------
static void runInteractive(const std::string &execPath) {
  // 1. Select a matrix market file and parse it into CSC.
  std::string filepath = SparseLib::UX::selectFileFromDirectory(execPath);
  std::cout << "\nProcessing (CSC) " << filepath << "...\n";
  Sparse::CSCMatrix dsm = SparseLib::IO::loadFromFile(filepath);

  // 2. Analyze original metrics.
  int fbm_before = SparseLib::Metrics::countFBM(dsm);
  long long tfbd_before = SparseLib::Metrics::fbmDiagonalDistance(dsm);
  std::cout << "FBM: " << fbm_before << ", TFBD: " << tfbd_before << "\n";

  // 3. Deduce SCCs using both algorithms.
  auto sccs_t = SparseLib::SCC::tarjanSCC(dsm);
  auto sccs_k = SparseLib::SCC::kosarajuSCC(dsm);
  std::cout << "SCCs (Tarjan): " << sccs_t.size()
            << ", (Kosaraju): " << sccs_k.size() << "\n";

  // 4. Build condensation graph.
  auto condense = SparseLib::SCC::condensationGraph(dsm, sccs_t);

  // 5. Topological sort on condensation graph.
  auto topo = SparseLib::SCC::topologicalSort(condense);

  // 6. Build permutation and apply P^T A P.
  auto perm = SparseLib::Permutation::buildPermutation(sccs_t, topo);
  auto permuted = SparseLib::Permutation::applyPermutation(dsm, perm);
  auto vr = SparseLib::Verify::verifyPermutation(dsm, permuted, perm);

  // Ensure the permutation is correct. i.e no information loss or introduced.
  if (!vr.passed)
    std::cerr << "Verification FAILED: " << vr.message << "\n";
  else
    std::cout << "Verification passed.\n";

  // Write permuted matrix into project out/permuted/ directory for Python
  // consumption
  try {
    fs::path base = fs::canonical(fs::path(execPath).parent_path() / "..");
    fs::path outDir = base / "out" / "permuted";
    fs::create_directories(outDir);
    fs::path inPath(filepath);
    fs::path outPath =
        outDir / (inPath.stem().string() + std::string("_permuted.mtx"));
    SparseLib::IO::saveToFile(permuted, outPath.string());
    std::cout << "Wrote permuted matrix to: " << outPath.string() << "\n";
  } catch (const std::exception &e) {
    std::cerr << "Warning: failed to write permuted matrix: " << e.what()
              << "\n";
  }

  // 7. Report before vs. after.
  int fbm_after = SparseLib::Metrics::countFBM(permuted);
  long long tfbd_after = SparseLib::Metrics::fbmDiagonalDistance(permuted);

  std::cout << "\nAfter BLTF permutation:\n"
            << "  FBM:  " << fbm_before << " -> " << fbm_after << "\n"
            << "  TFBD: " << tfbd_before << " -> " << tfbd_after << "\n";
}

// ---------------------------------------------------------------------------
// Benchmark mode: process all matrices, time every algorithm, export CSV
// ---------------------------------------------------------------------------
static void runBenchmark(const std::string &execPath) {
  SparseLib::Bench::printHeader();

  auto files = findMatrices(execPath);
  if (files.empty()) {
    std::cerr << "No .mtx files found in data/\n";
    return;
  }

  std::vector<SparseLib::Bench::BenchmarkRecord> records;

  for (const auto &filepath : files) {
    std::string name = fs::path(filepath).stem().string();
    Sparse::CSCMatrix dsm = SparseLib::IO::loadFromFile(filepath);

    // --- Time Tarjan SCC ---
    std::vector<std::vector<int>> sccs;
    auto t1 = SparseLib::Bench::measure(
        [&]() { sccs = SparseLib::SCC::tarjanSCC(dsm); });
    records.push_back({name, dsm.cols, dsm.nnz, "Tarjan", t1});

    // --- Time Kosaraju SCC ---
    std::vector<std::vector<int>> sccs_k;
    auto t2 = SparseLib::Bench::measure(
        [&]() { sccs_k = SparseLib::SCC::kosarajuSCC(dsm); });
    records.push_back({name, dsm.cols, dsm.nnz, "Kosaraju", t2});

    // --- Time condensation graph ---
    Sparse::CSCMatrix condense;
    auto t3 = SparseLib::Bench::measure(
        [&]() { condense = SparseLib::SCC::condensationGraph(dsm, sccs); });
    records.push_back({name, dsm.cols, dsm.nnz, "Condensation", t3});

    // --- Time topological sort ---
    std::vector<int> topo;
    auto t4 = SparseLib::Bench::measure(
        [&]() { topo = SparseLib::SCC::topologicalSort(condense); });
    records.push_back({name, dsm.cols, dsm.nnz, "TopoSort", t4});

    // --- Time permutation application ---
    auto perm = SparseLib::Permutation::buildPermutation(sccs, topo);
    Sparse::CSCMatrix permuted;
    auto t5 = SparseLib::Bench::measure([&]() {
      permuted = SparseLib::Permutation::applyPermutation(dsm, perm);
    });
    records.push_back({name, dsm.cols, dsm.nnz, "Permute", t5});

    SparseLib::Verify::Result vr;
    auto t6 = SparseLib::Bench::measure([&]() {
      vr = SparseLib::Verify::verifyPermutation(dsm, permuted, perm);
    });
    records.push_back({name, dsm.cols, dsm.nnz, "Verify", t6});
    if (!vr.passed) {
      std::cerr << "  Verification FAILED for " << name << ": " << vr.message
                << "\n";
      continue; // skip metrics for a broken permutation
    }
    std::cout << "  Verification assed.\n";

    // Write permuted matrix into project out/permuted/ directory (bench mode)
    try {
      fs::path base = fs::canonical(fs::path(execPath).parent_path() / "..");
      fs::path outDir = base / "out" / "permuted";
      fs::create_directories(outDir);
      fs::path inPath(filepath);
      fs::path outPath =
          outDir / (inPath.stem().string() + std::string("_permuted.mtx"));
      SparseLib::IO::saveToFile(permuted, outPath.string());
      std::cout << "  Wrote permuted matrix to: " << outPath.string() << "\n";
    } catch (const std::exception &e) {
      std::cerr << "  Warning: failed to write permuted matrix: " << e.what()
                << "\n";
    }

    // --- Quality metrics (before / after) ---
    int fbm_before = SparseLib::Metrics::countFBM(dsm);
    int fbm_after = SparseLib::Metrics::countFBM(permuted);
    long long tfbd_before = SparseLib::Metrics::fbmDiagonalDistance(dsm);
    long long tfbd_after = SparseLib::Metrics::fbmDiagonalDistance(permuted);

    int largest = 0, singletons = 0;
    for (const auto &c : sccs) {
      if (static_cast<int>(c.size()) > largest)
        largest = static_cast<int>(c.size());
      if (c.size() == 1)
        ++singletons;
    }

    SparseLib::Bench::MatrixStats stats;
    stats.name = name;
    stats.n = dsm.cols;
    stats.nnz = dsm.nnz;
    stats.fbm = fbm_before;
    stats.tfbd = tfbd_before;
    stats.num_sccs = static_cast<int>(sccs.size());
    stats.largest_scc = largest;
    stats.singleton_sccs = singletons;
    SparseLib::Bench::printMatrixStats(stats);

    std::cout << "    After:  FBM=" << fbm_after << ", TFBD=" << tfbd_after
              << "\n";
  }

  SparseLib::Bench::printTable(records);

  fs::path csvPath = fs::canonical(fs::path(execPath).parent_path() / "..") /
                     "out" / "benchmark.csv";
  fs::create_directories(csvPath.parent_path());
  SparseLib::Bench::writeCSV(csvPath.string(), records);
}

// ---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
  try {
    bool bench = false;
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]) == "--bench")
        bench = true;

    if (bench)
      runBenchmark(argv[0]);
    else
      runInteractive(argv[0]);

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
