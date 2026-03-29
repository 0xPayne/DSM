// CPSC 482 : DSM Project 
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/banding.hpp"
#include "../include/benchmark.hpp"
#include "../include/optimization_metrics.hpp"
#include "../include/parser.hpp"
#include "../include/permutation.hpp"
#include "../include/scc.hpp"
#include "../include/tearing.hpp"
#include "../include/ux.hpp"
#include "../include/verification.hpp"
#include "../include/output.hpp"

#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

static void reportQualityMetrics(const std::string &name,
                                 const Sparse::CSCMatrix &dsm,
                                 const std::vector<std::vector<int>> &sccs) {
  int fbm_before = SparseLib::Metrics::countFBM(dsm);
  long long tfbd_before = SparseLib::Metrics::fbmDiagonalDistance(dsm);

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
}

// Interactive mode
static void runInteractive(const std::string &execPath) {
  // Parse input
  std::string filepath = SparseLib::UX::selectFileFromDirectory(execPath);
  std::cout << "\nProcessing " << filepath << "...\n";
  Sparse::CSCMatrix dsm = SparseLib::IO::loadFromFile(filepath);

  // Baseline metrics
  int fbm_before = SparseLib::Metrics::countFBM(dsm);
  long long tfbd_before = SparseLib::Metrics::fbmDiagonalDistance(dsm);
  std::cout << "FBM: " << fbm_before << ", TFBD: " << tfbd_before << "\n";

  // SCC decomposition
  auto sccs = SparseLib::SCC::tarjanSCC(dsm);
  auto sccs_k = SparseLib::SCC::kosarajuSCC(dsm);
  std::cout << "SCCs (Tarjan): " << sccs.size()
            << ", (Kosaraju): " << sccs_k.size() << "\n";

  // Condensation and topological ordering
  auto condense = SparseLib::SCC::condensationGraph(dsm, sccs);
  auto topo = SparseLib::SCC::topologicalSort(condense);

  // Tearing and banding
  auto torn = SparseLib::Tearing::tearAllSCCs(sccs, dsm);
  auto banded = SparseLib::Banding::bandAllSCCs(torn, dsm);

  // Permutation
  auto perm = SparseLib::Permutation::buildPermutation(banded, topo);
  auto permuted = SparseLib::Permutation::applyPermutation(dsm, perm);

  // Verification
  auto vr = SparseLib::Verify::verifyPermutation(dsm, permuted, perm);
  if (!vr.passed)
    std::cerr << "Verification FAILED: " << vr.message << "\n";
  else
    std::cout << "Verification passed.\n";

  // Save result
  SparseLib::Output::savePermutedMatrix(execPath, filepath, permuted);
  SparseLib::Output::saveSCCBlocks(execPath, filepath, sccs, topo);

  // Report improvement
  int fbm_after = SparseLib::Metrics::countFBM(permuted);
  long long tfbd_after = SparseLib::Metrics::fbmDiagonalDistance(permuted);
  std::cout << "\nAfter BLTF permutation:\n"
            << "  FBM:  " << fbm_before << " -> " << fbm_after << "\n"
            << "  TFBD: " << tfbd_before << " -> " << tfbd_after << "\n";
}

// Benchmark mode: process all matrices, time every algorithm, export CSV
static void runBenchmark(const std::string &execPath) {
  SparseLib::Bench::printHeader();

  auto files = SparseLib::UX::listFiles(execPath);
  if (files.empty()) {
    std::cerr << "No .mtx files found in data/\n";
    return;
  }

  std::vector<SparseLib::Bench::BenchmarkRecord> records;
  std::vector<SparseLib::Bench::OptimizationRecord> optimization_records;

  for (const auto &filepath : files) {
    std::string name = fs::path(filepath).stem().string();
    Sparse::CSCMatrix dsm = SparseLib::IO::loadFromFile(filepath);

    // Time each thing
    std::vector<std::vector<int>> sccs;
    auto t1 = SparseLib::Bench::measure(
        [&]() { sccs = SparseLib::SCC::tarjanSCC(dsm); });
    records.push_back({name, dsm.cols, dsm.nnz, "Tarjan", t1});

    std::vector<std::vector<int>> sccs_k;
    auto t2 = SparseLib::Bench::measure(
        [&]() { sccs_k = SparseLib::SCC::kosarajuSCC(dsm); });
    records.push_back({name, dsm.cols, dsm.nnz, "Kosaraju", t2});

    Sparse::CSCMatrix condense;
    auto t3 = SparseLib::Bench::measure(
        [&]() { condense = SparseLib::SCC::condensationGraph(dsm, sccs); });
    records.push_back({name, dsm.cols, dsm.nnz, "Condensation", t3});

    std::vector<int> topo;
    auto t4 = SparseLib::Bench::measure(
        [&]() { topo = SparseLib::SCC::topologicalSort(condense); });
    records.push_back({name, dsm.cols, dsm.nnz, "TopoSort", t4});

    std::vector<std::vector<int>> torn_sccs;
    auto t_tear = SparseLib::Bench::measure(
        [&]() { torn_sccs = SparseLib::Tearing::tearAllSCCs(sccs, dsm); });
    records.push_back({name, dsm.cols, dsm.nnz, "Tearing", t_tear});

    std::vector<std::vector<int>> banded_sccs;
    auto t_band = SparseLib::Bench::measure([&]() {
      banded_sccs = SparseLib::Banding::bandAllSCCs(torn_sccs, dsm);
    });
    records.push_back({name, dsm.cols, dsm.nnz, "Banding", t_band});

    // Recompute stage orderings once for quality metrics export.
    torn_sccs = SparseLib::Tearing::tearAllSCCs(sccs, dsm);
    banded_sccs = SparseLib::Banding::bandAllSCCs(torn_sccs, dsm);

    auto perm = SparseLib::Permutation::buildPermutation(banded_sccs, topo);
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
      continue;
    }
    std::cout << "  Verification Passed.\n";

    // Save and report
    SparseLib::Output::savePermutedMatrix(execPath, filepath, permuted);
    SparseLib::Output::saveSCCBlocks(execPath, filepath, sccs, topo);
    reportQualityMetrics(name, dsm, sccs);

    auto perm_topo = SparseLib::Permutation::buildPermutation(sccs, topo);
    auto perm_tearing = SparseLib::Permutation::buildPermutation(torn_sccs, topo);
    auto topo_only = SparseLib::Permutation::applyPermutation(dsm, perm_topo);
    auto torn_only = SparseLib::Permutation::applyPermutation(dsm, perm_tearing);

    SparseLib::Bench::OptimizationRecord quality;
    quality.matrix_name = name;
    quality.n = dsm.cols;
    quality.nnz = dsm.nnz;
    quality.fbm_base = SparseLib::Metrics::countFBM(dsm);
    quality.fbm_topo = SparseLib::Metrics::countFBM(topo_only);
    quality.fbm_tearing = SparseLib::Metrics::countFBM(torn_only);
    quality.fbm_banding = SparseLib::Metrics::countFBM(permuted);
    quality.tfbd_base = SparseLib::Metrics::fbmDiagonalDistance(dsm);
    quality.tfbd_topo = SparseLib::Metrics::fbmDiagonalDistance(topo_only);
    quality.tfbd_tearing = SparseLib::Metrics::fbmDiagonalDistance(torn_only);
    quality.tfbd_banding = SparseLib::Metrics::fbmDiagonalDistance(permuted);
    optimization_records.push_back(quality);
    SparseLib::Bench::printOptimizationComparison(quality);
  }

  SparseLib::Bench::printTable(records);

  fs::path csvPath =
      SparseLib::UX::projectRoot(execPath) / "out" / "benchmark.csv";
  fs::create_directories(csvPath.parent_path());
  SparseLib::Bench::writeCSV(csvPath.string(), records);

  fs::path optimizationPath =
      SparseLib::UX::projectRoot(execPath) / "out" / "optimization_metrics.csv";
  SparseLib::Bench::writeOptimizationCSV(optimizationPath.string(), optimization_records);
}

// Entry point
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
