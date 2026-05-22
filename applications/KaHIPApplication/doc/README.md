# KaHIP Documentation

**KaHIP — Karlsruhe High Quality Partitioning**
Version 3.25 · MIT License · C++11/14

This documentation covers the KaHIP graph partitioning framework in depth: its architecture, algorithms, APIs, build system, configuration options, and usage patterns.

---

## Document Index

| Document                                 | Description                                                  |
|------------------------------------------|--------------------------------------------------------------|
| [architecture.md](architecture.md)       | Library architecture, module layout, data flow               |
| [build.md](build.md)                     | Building from source: CMake options, dependencies, platforms |
| [api_serial.md](api_serial.md)           | Serial C/C++ API reference (kaHIP_interface.h)               |
| [api_parallel.md](api_parallel.md)       | Parallel API reference (parhip_interface.h, kaffpaE)         |
| [algorithms.md](algorithms.md)           | Algorithm descriptions: KaFFPa, ParHIP, KaBaPE, etc.         |
| [graph_format.md](graph_format.md)       | Input/output graph file formats (METIS, BGF)                 |
| [configuration.md](configuration.md)     | All configuration parameters and preconfigurations           |
| [python_bindings.md](python_bindings.md) | Python interface and usage                                   |
| [examples.md](examples.md)               | Complete usage examples for all major features               |

---

## Quick Start

```bash
# Build (from repository root)
cd build && ./configure.sh

# Partition a graph into 4 parts (ECO quality)
./deploy/kaffpa examples/rgg_n_2_15_s0.graph --k 4 --preconfiguration=eco

# Use the library from C++
#include "kaHIP_interface.h"
kaffpa(&n, vwgt, xadj, adjcwgt, adjncy, &nparts, &imbalance, false, 0, ECO, &edgecut, part);

# Parallel partitioning
mpirun -n 8 ./deploy/parhip examples/rgg_n_2_15_s0.graph --k 32 --preconfiguration=fastmesh
```

---

## What is KaHIP?

KaHIP solves the **graph partitioning problem**: given a graph G = (V, E) with node weights w: V → ℝ⁺ and edge weights c: E → ℝ⁺, find a partition of V into k equally-weighted blocks V₁, …, Vₖ that minimizes the total weight of edges crossing block boundaries (the **edge cut**), subject to a balance constraint:

```
w(Vᵢ) ≤ (1 + ε) · W/k    for all i ∈ {1, …, k}
```

where W = Σ w(v) is the total node weight and ε is the allowed imbalance.

This problem is NP-hard in general. KaHIP uses multilevel heuristics that achieve near-optimal quality in practice.

### Applications

- **Parallel scientific computing**: domain decomposition for FEM/FVM solvers
- **Distributed graph processing**: workload balancing across compute nodes
- **VLSI design**: circuit partitioning
- **Social network analysis**: community detection
- **Sparse matrix reordering**: fill-in minimization for direct solvers

### Algorithm Family

| Component   | Description                                                        |
|-------------|--------------------------------------------------------------------|
| **KaFFPa**  | Karlsruhe Fast Flow Partitioner — multilevel serial partitioner    |
| **KaFFPaE** | KaFFPa Evolutionary — MPI-parallel memetic algorithm               |
| **KaBaPE**  | Extension of KaFFPaE with additional combine operators             |
| **ParHIP**  | Parallel Hierarchical Partitioner — distributed-memory partitioner |
| **Buffoon** | Specialized for road networks                                      |

---

## Repository Layout

```
KaHIP/
├── app/                    # Main application entry points (executables)
├── build/                  # Build output directory
│   └── configure.sh        # Comprehensive build script (this repo)
├── cmake/                  # CMake modules (FindGurobi.cmake)
├── doc/                    # This documentation
├── examples/               # Example graph files
├── extern/                 # Bundled dependencies (argtable3)
├── interface/              # Public C/C++ serial API
│   ├── kaHIP_interface.h   # Serial header
│   └── kaHIP_interface.cpp # Serial implementation
├── lib/                    # Core library implementation
│   ├── algorithms/         # Graph algorithms (SCC, push-relabel, …)
│   ├── data_structure/     # Graph, priority queues, matrices
│   ├── io/                 # Graph I/O
│   ├── mapping/            # Process mapping
│   ├── node_ordering/      # Nested dissection ordering
│   ├── parallel_mh/        # Parallel memetic algorithm
│   ├── partition/          # Partitioning: coarsening, IP, refinement
│   ├── spac/               # Edge partitioning (SPAC)
│   └── tools/              # Utilities, metrics, random functions
├── misc/                   # Examples, Python module, Java bindings
│   ├── example_library_call/
│   ├── example_parhip_call/
│   ├── java_jni_wrapper/
│   └── pymodule/
├── parallel/               # MPI/parallel implementations
│   ├── modified_kahip/     # KaHIP variant used by ParHIP internally
│   └── parallel_src/       # ParHIP distributed partitioner
│       └── interface/      # parhip_interface.h
├── python/                 # Python package source
│   └── kahip/
│       ├── __init__.py
│       └── _version.py
├── CMakeLists.txt          # Root CMake build file
└── README.md               # Project README
```
