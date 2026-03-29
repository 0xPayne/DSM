// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once

#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace SparseLib::Bench {

static constexpr int DEFAULT_WARMUP = 2;
static constexpr int DEFAULT_TRIALS = 5;

struct TimingResult {
    std::vector<double> samples_us;
    double min_us    = 0.0;
    double max_us    = 0.0;
    double mean_us   = 0.0;
    double median_us = 0.0;
    double stddev_us = 0.0;
};

// Compiler optimization barrier: prevents the compiler from reordering
// memory operations across this point. Same technique as Google Benchmark
// and folly::doNotOptimizeAway. Ensures timing measurements aren't
// polluted by instruction reordering.
#if defined(__GNUC__) || defined(__clang__)
inline void clobberMemory() { asm volatile("" : : : "memory"); }
#else
inline void clobberMemory() { std::atomic_signal_fence(std::memory_order_seq_cst); }
#endif

// Benchmark any callable over trials timed runs preceded by warmup
// untimed runs.
template<typename Func>
TimingResult measure(Func&& fn, int warmup = DEFAULT_WARMUP, int trials = DEFAULT_TRIALS) {
    for (int i = 0; i < warmup; ++i) fn();

    std::vector<double> samples(trials);
    for (int i = 0; i < trials; ++i) {
        clobberMemory();
        auto t0 = std::chrono::steady_clock::now();
        fn();
        clobberMemory();
        auto t1 = std::chrono::steady_clock::now();
        samples[i] = std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    auto sorted = samples;
    std::sort(sorted.begin(), sorted.end());

    double sum  = std::accumulate(sorted.begin(), sorted.end(), 0.0);
    double mean = sum / trials;

    double var = 0.0;
    for (double s : sorted) var += (s - mean) * (s - mean);
    double stddev = (trials > 1) ? std::sqrt(var / (trials - 1)) : 0.0;

    double median = (trials % 2 == 0)
        ? (sorted[trials / 2 - 1] + sorted[trials / 2]) / 2.0
        : sorted[trials / 2];

    return {samples, sorted.front(), sorted.back(), mean, median, stddev};
}

struct BenchmarkRecord {
    std::string matrix_name;
    int n   = 0;
    int nnz = 0;
    std::string algorithm;
    TimingResult timing;
};

struct MatrixStats {
    std::string name;
    int n          = 0;
    int nnz        = 0;
    double density = 0.0;
    int fbm        = 0;
    long long tfbd = 0;
    int num_sccs   = 0;
    int largest_scc = 0;
    int singleton_sccs = 0;
};

// Per-matrix optimization comparison across reordering stages.
struct OptimizationRecord {
    std::string matrix_name;
    int n   = 0;
    int nnz = 0;
    int fbm_base    = 0;
    int fbm_topo    = 0;
    int fbm_tearing = 0;
    int fbm_banding = 0;
    long long tfbd_base    = 0;
    long long tfbd_topo    = 0;
    long long tfbd_tearing = 0;
    long long tfbd_banding = 0;
};

std::string systemInfo();
void printHeader();
void printMatrixStats(const MatrixStats& stats);
void printOptimizationComparison(const OptimizationRecord& record);
void printTable(const std::vector<BenchmarkRecord>& records);
void writeCSV(const std::string& path, const std::vector<BenchmarkRecord>& records);
void writeOptimizationCSV(const std::string& path, const std::vector<OptimizationRecord>& records);

} // namespace SparseLib::Bench
