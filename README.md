# DSM
CPSC 482 - Design Strucuture Matricies Optimization Project.


## Project Outline
```text
.
├── Makefile
├── LICENSE
├── README.md
├── data/            # Sample data (Benchmark DSMs)
├── include/         # Header files (.hpp)
├── lib/             # External libraries
├── out/             # Compiled binaries
├── report/          # LaTeX source files
└── src/             # Source code (.cpp)
```

## Run Benchmarks and Export Report Metrics
```bash
make
./out/main --bench
```

Generated outputs:
- `out/benchmark.csv` (algorithm timing metrics)
- `out/optimization_metrics.csv` (baseline vs topo-only vs tearing vs banding quality metrics)
- `out/permuted/*_permuted.mtx` (final permuted matrices)

## Plot Numerical Testing Results
```bash
python3 scripts/plot_benchmarks.py
python3 scripts/plot_optimization_metrics.py
python3 scripts/generate_report_tables.py
```
