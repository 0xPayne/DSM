// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/benchmark.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>

#ifdef __APPLE__
#include <sys/sysctl.h>
#include <sys/types.h>
#elif defined(__linux__)
#include <fstream>
#endif

namespace SparseLib::Bench {

// System information: CPU, memory, cache hierarchy, compiler

#ifdef __APPLE__
static std::string sysctlString(const char* key) {
    char buf[256] = {};
    size_t len = sizeof(buf);
    if (sysctlbyname(key, buf, &len, nullptr, 0) == 0) return buf;
    return "";
}

static int64_t sysctlInt64(const char* key) {
    int64_t val = 0;
    size_t len = sizeof(val);
    if (sysctlbyname(key, &val, &len, nullptr, 0) == 0) return val;
    return -1;
}
#endif

std::string systemInfo() {
    std::ostringstream oss;

#ifdef __APPLE__
    std::string brand = sysctlString("machdep.cpu.brand_string");
    if (brand.empty()) brand = sysctlString("hw.machine");
    if (!brand.empty()) oss << "  CPU:       " << brand << "\n";

    int64_t mem = sysctlInt64("hw.memsize");
    if (mem > 0) oss << "  RAM:       " << (mem / (1024 * 1024 * 1024)) << " GB\n";

    // Intel keys first, then Apple Silicon fallback
    int64_t l1d = sysctlInt64("hw.l1dcachesize");
    if (l1d <= 0) l1d = sysctlInt64("hw.perflevel0.l1dcachesize");
    if (l1d > 0) oss << "  L1d:       " << (l1d / 1024) << " KB\n";

    int64_t l1i = sysctlInt64("hw.l1icachesize");
    if (l1i <= 0) l1i = sysctlInt64("hw.perflevel0.l1icachesize");
    if (l1i > 0) oss << "  L1i:       " << (l1i / 1024) << " KB\n";

    int64_t l2 = sysctlInt64("hw.l2cachesize");
    if (l2 <= 0) l2 = sysctlInt64("hw.perflevel0.l2cachesize");
    if (l2 > 0) {
        if (l2 >= 1024 * 1024)
            oss << "  L2:        " << (l2 / (1024 * 1024)) << " MB\n";
        else
            oss << "  L2:        " << (l2 / 1024) << " KB\n";
    }

    int64_t l3 = sysctlInt64("hw.l3cachesize");
    if (l3 > 0) oss << "  L3:        " << (l3 / (1024 * 1024)) << " MB\n";

    int64_t cacheline = sysctlInt64("hw.cachelinesize");
    if (cacheline > 0) oss << "  Cacheline: " << cacheline << " B\n";

    int64_t pcores = sysctlInt64("hw.perflevel0.physicalcpu");
    int64_t ecores = sysctlInt64("hw.perflevel1.physicalcpu");
    if (pcores > 0) oss << "  P-cores:   " << pcores << "\n";
    if (ecores > 0) oss << "  E-cores:   " << ecores << "\n";

#elif defined(__linux__)
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line)) {
        if (line.find("model name") != std::string::npos) {
            oss << "  CPU:       " << line.substr(line.find(':') + 2) << "\n";
            break;
        }
    }

    std::ifstream meminfo("/proc/meminfo");
    while (std::getline(meminfo, line)) {
        if (line.find("MemTotal") != std::string::npos) {
            oss << "  RAM:       " << line.substr(line.find(':') + 2) << "\n";
            break;
        }
    }

    for (int i = 0; i <= 3; ++i) {
        std::string base = "/sys/devices/system/cpu/cpu0/cache/index" + std::to_string(i);
        std::ifstream sz(base + "/size"), tp(base + "/type");
        if (sz && tp) {
            std::string size_str, type_str;
            std::getline(tp, type_str);
            std::getline(sz, size_str);
            oss << "  L" << i << " (" << type_str << "): " << size_str << "\n";
        }
    }
#endif

    oss << "  Compiler:  " << __VERSION__ << "\n";
    oss << "  Flags:     "
#ifdef __OPTIMIZE__
        "optimized (-O)"
#else
        "unoptimized"
#endif
        "\n";

    return oss.str();
}

// Output formatting

static const char* SEPARATOR =
    "================================================================================";

static double reductionPercent(long long base, long long stage) {
    if (base <= 0) return 0.0;
    return 100.0 * static_cast<double>(base - stage) / static_cast<double>(base);
}

void printHeader() {
    std::cout << "\n" << SEPARATOR << "\n"
              << "  DSM Benchmark Suite\n"
              << SEPARATOR << "\n\n"
              << "System:\n" << systemInfo() << "\n"
              << "Parameters:\n"
              << "  Warm-up runs: " << DEFAULT_WARMUP << "\n"
              << "  Timed runs:   " << DEFAULT_TRIALS << "\n\n"
              << SEPARATOR << "\n";
}

void printMatrixStats(const MatrixStats& stats) {
    double pct = (stats.n > 0)
        ? 100.0 * stats.nnz / ((double)stats.n * stats.n)
        : 0.0;

    std::cout << "\n  " << stats.name
              << "  (n=" << stats.n
              << ", nnz=" << stats.nnz
              << ", density=" << std::fixed << std::setprecision(4) << pct << "%)\n"
              << "    FBM: " << stats.fbm
              << "  |  TFBD: " << stats.tfbd
              << "  |  SCCs: " << stats.num_sccs
              << " [largest: " << stats.largest_scc
              << ", singletons: " << stats.singleton_sccs << "]\n";
}

void printOptimizationComparison(const OptimizationRecord& record) {
    std::cout << "    Stage quality (FBM / TFBD):\n"
              << "      Baseline : " << record.fbm_base << " / " << record.tfbd_base << "\n"
              << "      TopoOnly : " << record.fbm_topo << " / " << record.tfbd_topo
              << "  (" << std::fixed << std::setprecision(2)
              << reductionPercent(record.fbm_base, record.fbm_topo) << "% FBM, "
              << reductionPercent(record.tfbd_base, record.tfbd_topo) << "% TFBD)\n"
              << "      Tearing  : " << record.fbm_tearing << " / " << record.tfbd_tearing
              << "  (" << reductionPercent(record.fbm_base, record.fbm_tearing) << "% FBM, "
              << reductionPercent(record.tfbd_base, record.tfbd_tearing) << "% TFBD)\n"
              << "      Banding  : " << record.fbm_banding << " / " << record.tfbd_banding
              << "  (" << reductionPercent(record.fbm_base, record.fbm_banding) << "% FBM, "
              << reductionPercent(record.tfbd_base, record.tfbd_banding) << "% TFBD)\n";
}

void printTable(const std::vector<BenchmarkRecord>& records) {
    if (records.empty()) return;

    std::cout << "\n" << SEPARATOR << "\n"
              << "  Timing Results\n"
              << SEPARATOR << "\n\n";

    // Column widths
    const int WN = 20, WI = 8, WZ = 10, WA = 14, WS = 12, WR = 10;

    // Header
    std::printf("%-*s %*s %*s %-*s %*s %*s %*s %*s %*s %*s %*s\n",
        WN, "Matrix", WI, "n", WZ, "nnz", WA, "Algorithm",
        WS, "Mean(us)", WS, "Std(us)", WS, "Min(us)",
        WS, "Median(us)", WS, "Max(us)", WR, "ns/vtx", WR, "ns/edge");

    int total_w = WN + WI + WZ + WA + 5*WS + 2*WR + 10;
    std::cout << std::string(total_w, '-') << "\n";

    for (const auto& r : records) {
        double ns_vtx  = (r.n > 0)   ? (r.timing.mean_us * 1000.0 / r.n)   : 0.0;
        double ns_edge = (r.nnz > 0) ? (r.timing.mean_us * 1000.0 / r.nnz) : 0.0;

        std::printf("%-*s %*d %*d %-*s %*.2f %*.2f %*.2f %*.2f %*.2f %*.2f %*.2f\n",
            WN, r.matrix_name.c_str(),
            WI, r.n,
            WZ, r.nnz,
            WA, r.algorithm.c_str(),
            WS, r.timing.mean_us,
            WS, r.timing.stddev_us,
            WS, r.timing.min_us,
            WS, r.timing.median_us,
            WS, r.timing.max_us,
            WR, ns_vtx,
            WR, ns_edge);
    }

    std::cout << std::string(total_w, '-') << "\n";
}

// CSV export

void writeCSV(const std::string& path, const std::vector<BenchmarkRecord>& records) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to open " << path << " for writing.\n";
        return;
    }

    int num_runs = DEFAULT_TRIALS;
    if (!records.empty())
        num_runs = static_cast<int>(records[0].timing.samples_us.size());

    out << "matrix,n,nnz,algorithm,mean_us,stddev_us,min_us,median_us,max_us,"
        << "ns_per_vertex,ns_per_edge";
    for (int i = 0; i < num_runs; ++i) out << ",run" << (i + 1) << "_us";
    out << "\n";

    out << std::fixed << std::setprecision(2);
    for (const auto& r : records) {
        double ns_vtx  = (r.n > 0)   ? (r.timing.mean_us * 1000.0 / r.n)   : 0.0;
        double ns_edge = (r.nnz > 0) ? (r.timing.mean_us * 1000.0 / r.nnz) : 0.0;

        out << r.matrix_name << ","
            << r.n << "," << r.nnz << ","
            << r.algorithm << ","
            << r.timing.mean_us << ","
            << r.timing.stddev_us << ","
            << r.timing.min_us << ","
            << r.timing.median_us << ","
            << r.timing.max_us << ","
            << ns_vtx << "," << ns_edge;

        for (double s : r.timing.samples_us) out << "," << s;
        out << "\n";
    }

    std::cout << "\nCSV written to: " << path << "\n";
}

void writeOptimizationCSV(const std::string& path, const std::vector<OptimizationRecord>& records) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to open " << path << " for writing.\n";
        return;
    }

    out << "matrix,n,nnz,density_pct,"
        << "fbm_base,fbm_topo,fbm_tearing,fbm_banding,"
        << "tfbd_base,tfbd_topo,tfbd_tearing,tfbd_banding,"
        << "fbm_gain_topo_abs,fbm_gain_tearing_abs,fbm_gain_banding_abs,"
        << "tfbd_gain_topo_abs,tfbd_gain_tearing_abs,tfbd_gain_banding_abs,"
        << "fbm_gain_topo_pct,fbm_gain_tearing_pct,fbm_gain_banding_pct,"
        << "tfbd_gain_topo_pct,tfbd_gain_tearing_pct,tfbd_gain_banding_pct,"
        << "fbm_gain_topo_to_tearing_abs,fbm_gain_tearing_to_banding_abs,"
        << "tfbd_gain_topo_to_tearing_abs,tfbd_gain_tearing_to_banding_abs\n";

    out << std::fixed << std::setprecision(4);
    for (const auto& r : records) {
        double density_pct = (r.n > 0)
            ? 100.0 * static_cast<double>(r.nnz) / (static_cast<double>(r.n) * r.n)
            : 0.0;

        int fbm_gain_topo = r.fbm_base - r.fbm_topo;
        int fbm_gain_tear = r.fbm_base - r.fbm_tearing;
        int fbm_gain_band = r.fbm_base - r.fbm_banding;
        long long tfbd_gain_topo = r.tfbd_base - r.tfbd_topo;
        long long tfbd_gain_tear = r.tfbd_base - r.tfbd_tearing;
        long long tfbd_gain_band = r.tfbd_base - r.tfbd_banding;

        out << r.matrix_name << ","
            << r.n << "," << r.nnz << "," << density_pct << ","
            << r.fbm_base << "," << r.fbm_topo << "," << r.fbm_tearing << "," << r.fbm_banding << ","
            << r.tfbd_base << "," << r.tfbd_topo << "," << r.tfbd_tearing << "," << r.tfbd_banding << ","
            << fbm_gain_topo << "," << fbm_gain_tear << "," << fbm_gain_band << ","
            << tfbd_gain_topo << "," << tfbd_gain_tear << "," << tfbd_gain_band << ","
            << reductionPercent(r.fbm_base, r.fbm_topo) << ","
            << reductionPercent(r.fbm_base, r.fbm_tearing) << ","
            << reductionPercent(r.fbm_base, r.fbm_banding) << ","
            << reductionPercent(r.tfbd_base, r.tfbd_topo) << ","
            << reductionPercent(r.tfbd_base, r.tfbd_tearing) << ","
            << reductionPercent(r.tfbd_base, r.tfbd_banding) << ","
            << (r.fbm_topo - r.fbm_tearing) << "," << (r.fbm_tearing - r.fbm_banding) << ","
            << (r.tfbd_topo - r.tfbd_tearing) << "," << (r.tfbd_tearing - r.tfbd_banding)
            << "\n";
    }

    std::cout << "Optimization CSV written to: " << path << "\n";
}

void writeOptimizationSummaryCSV(const std::string& path, const std::vector<OptimizationRecord>& records) {
    std::vector<OptimizationRecord> filtered;
    filtered.reserve(records.size());
    for (const auto& r : records) {
        if (r.matrix_name != "easy-example")
            filtered.push_back(r);
    }
    if (filtered.empty()) {
        std::cerr << "No optimization records for summary (after excluding easy-example).\n";
        return;
    }

    double sum_fbm_topo_pct = 0, sum_fbm_tear_pct = 0, sum_fbm_band_pct = 0;
    double sum_tfbd_topo_pct = 0, sum_tfbd_tear_pct = 0, sum_tfbd_band_pct = 0;
    long long sum_fbm_base = 0, sum_fbm_topo = 0, sum_fbm_tear = 0, sum_fbm_band = 0;
    long long sum_tfbd_base = 0, sum_tfbd_topo = 0, sum_tfbd_tear = 0, sum_tfbd_band = 0;

    for (const auto& r : filtered) {
        sum_fbm_topo_pct += reductionPercent(r.fbm_base, r.fbm_topo);
        sum_fbm_tear_pct += reductionPercent(r.fbm_base, r.fbm_tearing);
        sum_fbm_band_pct += reductionPercent(r.fbm_base, r.fbm_banding);
        sum_tfbd_topo_pct += reductionPercent(r.tfbd_base, r.tfbd_topo);
        sum_tfbd_tear_pct += reductionPercent(r.tfbd_base, r.tfbd_tearing);
        sum_tfbd_band_pct += reductionPercent(r.tfbd_base, r.tfbd_banding);

        sum_fbm_base += r.fbm_base;
        sum_fbm_topo += r.fbm_topo;
        sum_fbm_tear += r.fbm_tearing;
        sum_fbm_band += r.fbm_banding;
        sum_tfbd_base += r.tfbd_base;
        sum_tfbd_topo += r.tfbd_topo;
        sum_tfbd_tear += r.tfbd_tearing;
        sum_tfbd_band += r.tfbd_banding;
    }

    const int n = static_cast<int>(filtered.size());
    const double inv_n = 1.0 / static_cast<double>(n);

    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to open " << path << " for writing.\n";
        return;
    }

    out << std::fixed << std::setprecision(4);
    out << "metric,value\n";
    out << "matrices_included," << n << "\n";
    out << "mean_fbm_gain_topo_pct," << (sum_fbm_topo_pct * inv_n) << "\n";
    out << "mean_fbm_gain_tearing_pct," << (sum_fbm_tear_pct * inv_n) << "\n";
    out << "mean_fbm_gain_banding_pct," << (sum_fbm_band_pct * inv_n) << "\n";
    out << "mean_tfbd_gain_topo_pct," << (sum_tfbd_topo_pct * inv_n) << "\n";
    out << "mean_tfbd_gain_tearing_pct," << (sum_tfbd_tear_pct * inv_n) << "\n";
    out << "mean_tfbd_gain_banding_pct," << (sum_tfbd_band_pct * inv_n) << "\n";
    out << "macro_fbm_gain_topo_pct," << reductionPercent(sum_fbm_base, sum_fbm_topo) << "\n";
    out << "macro_fbm_gain_tearing_pct," << reductionPercent(sum_fbm_base, sum_fbm_tear) << "\n";
    out << "macro_fbm_gain_banding_pct," << reductionPercent(sum_fbm_base, sum_fbm_band) << "\n";
    out << "macro_tfbd_gain_topo_pct," << reductionPercent(sum_tfbd_base, sum_tfbd_topo) << "\n";
    out << "macro_tfbd_gain_tearing_pct," << reductionPercent(sum_tfbd_base, sum_tfbd_tear) << "\n";
    out << "macro_tfbd_gain_banding_pct," << reductionPercent(sum_tfbd_base, sum_tfbd_band) << "\n";

    std::cout << "Optimization summary CSV written to: " << path << "\n";
    std::cout << "  (mean over " << n << " matrices, excluding easy-example)\n";
    std::cout << "  Mean FBM%% vs baseline — topo / tear / band: "
              << (sum_fbm_topo_pct * inv_n) << " / " << (sum_fbm_tear_pct * inv_n) << " / "
              << (sum_fbm_band_pct * inv_n) << "\n";
    std::cout << "  Mean TFBD%% vs baseline — topo / tear / band: "
              << (sum_tfbd_topo_pct * inv_n) << " / " << (sum_tfbd_tear_pct * inv_n) << " / "
              << (sum_tfbd_band_pct * inv_n) << "\n";
}

} // namespace SparseLib::Bench
