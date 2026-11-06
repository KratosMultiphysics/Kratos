# KaHIP Architecture

## High-Level Architecture

KaHIP is organized as a layered library with public C interfaces wrapping a rich internal C++ implementation.

```
┌─────────────────────────────────────────────────────────────────┐
│                        User Applications                        │
├────────────────────┬────────────────────┬───────────────────────┤
│   C/C++ Interface  │  Python Interface  │  CLI Executables      │
│  kaHIP_interface.h │  kahip (PyPI)      │  kaffpa, parhip, …    │
├────────────────────┴────────────────────┴───────────────────────┤
│                    Partitioning Engine                          │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │                   graph_partitioner                       │  │
│  │  ┌──────────────┐  ┌────────────────┐  ┌──────────────┐   │  │
│  │  │  Coarsening  │→ │    Initial     │→ │ Uncoarsening │   │  │
│  │  │  Phase       │  │  Partitioning  │  │ + Refinement │   │  │
│  │  └──────────────┘  └────────────────┘  └──────────────┘   │  │
│  └───────────────────────────────────────────────────────────┘  │
├─────────────────────────────────────────────────────────────────┤
│              Supporting Library Modules                         │
│   algorithms/   data_structure/   tools/   io/   mapping/       │
│   node_ordering/   spac/   parallel_mh/                         │
├─────────────────────────────────────────────────────────────────┤
│            Parallel Layer (MPI)                                 │
│   parallel_mh/   parallel/parallel_src/   lib/tools/mpi_tools   │
└─────────────────────────────────────────────────────────────────┘
```

---

## Module Breakdown

### `interface/`

The public-facing C API. This is the only module users need to link against.

- **`kaHIP_interface.h`** — declares all public functions and mode constants
- **`kaHIP_interface.cpp`** — translates the flat C interface into internal C++ calls, sets up `PartitionConfig`, invokes `graph_partitioner`

All graph data structures exposed in the interface use CSR (Compressed Sparse Row) format, matching the METIS API conventions.

### `lib/partition/`

The heart of KaHIP. Implements the multilevel partitioning framework.

```
lib/partition/
├── partition_config.h          # 100+ tunable parameters (struct PartitionConfig)
├── graph_partitioner.{h,cpp}   # Top-level orchestrator
├── w_cycles/
│   └── wcycle_partitioner      # W-cycle vs V-cycle multilevel scheme
├── coarsening/                 # Phase 1: graph compression
│   ├── coarsening.{h,cpp}      # Main coarsening loop
│   ├── contraction.{h,cpp}     # Node merging
│   ├── edge_rating/            # Edge importance metrics
│   │   └── edge_ratings        # EXPANSIONSTAR, WEIGHT, PSEUDOGEOM, …
│   ├── matching/               # Matching algorithms
│   │   ├── random_matching     # Random matching
│   │   ├── gpa/                # Global Path Algorithm (optimal matching)
│   │   └── gpa_matching        # RANDOM_GPA (randomized GPA)
│   └── clustering/
│       └── size_constraint_label_propagation  # CLUSTER_COARSENING
├── initial_partitioning/       # Phase 2: partition coarsest graph
│   ├── initial_partitioning    # Dispatcher
│   ├── initial_partitioner     # Tries multiple methods, keeps best
│   ├── bipartition             # BFS-based 2-way split
│   └── initial_refinement/     # FM refinement at coarsest level
└── uncoarsening/               # Phase 3: project and refine
    ├── uncoarsening.{h,cpp}    # Project partition up hierarchy
    ├── separator/              # Node separator computation
    │   ├── vertex_separator_algorithm
    │   └── vertex_separator_flow_solver
    └── refinement/             # Local search improvements
        ├── refinement.{h,cpp}  # Dispatcher
        ├── mixed_refinement    # 2-way + k-way
        ├── label_propagation_refinement/
        ├── quotient_graph_refinement/   # 2-way FM + flow
        │   ├── 2way_fm_refinement/      # Fiduccia-Mattheyses
        │   └── flow_refinement/         # Max-flow based
        ├── kway_graph_refinement/       # k-way FM variants
        ├── cycle_improvements/          # KaBaPE cycle refinement
        ├── node_separators/             # NS-specific refinement
        └── tabu_search/                 # Tabu search
```

### `lib/data_structure/`

Core graph representation used throughout the library.

- **`graph_access`** — the main graph class; stores adjacency in CSR, node/edge weights, partition labels
- **`graph_hierarchy`** — stack of coarsened graphs built during coarsening, used for projection
- **`priority_queues/`** — bucket queues and binary heaps used by FM and flow refinement
- **`matrix/`** — sparse matrix types for process mapping QAP

### `lib/algorithms/`

Standalone graph algorithms:

| File | Algorithm |
|---|---|
| `strongly_connected_components` | Kosaraju/Tarjan SCC |
| `topological_sort` | Topological ordering |
| `push_relabel` | Max-flow via push-relabel (used by flow refinement) |
| `cycle_search` | Negative cycle detection (used by KaBaPE) |

### `lib/tools/`

Utility functions used across modules:

| File | Purpose |
|---|---|
| `quality_metrics` | Compute edge cut, imbalance, connectivity |
| `random_functions` | Seeded PRNG, permutation generation |
| `graph_extractor` | Extract subgraphs and quotient graphs |
| `partition_snapshooter` | Save/load intermediate partitions |
| `graph_communication` | MPI graph communication primitives |
| `mpi_tools` | MPI utility wrappers |
| `misc` | Miscellaneous helpers |

### `lib/io/`

Graph file I/O:

- **`graph_io`** — reads/writes METIS format (`.graph`), binary format (`.bgf`), and partition files
- Supports weighted and unweighted graphs
- Optional `--mmap_io` for large files

### `lib/mapping/`

Process mapping algorithms for minimizing communication time on hierarchical processor networks. Models the problem as a Quadratic Assignment Problem (QAP).

| File | Description |
|---|---|
| `local_search_mapping` | Local search for QAP |
| `fast_construct_mapping` | Fast initial mapping construction |
| `mapping_algorithms` | Top-level dispatcher |
| `construct_distance_matrix` | Build processor distance matrix |
| `full_search_space` | Enumerate all mapping assignments |
| `communication_graph_search_space` | Restrict search to communication graph |

### `lib/node_ordering/`

Node ordering algorithms for fill-in minimization in sparse direct solvers.

| File | Description |
|---|---|
| `nested_dissection` | Multilevel nested dissection |
| `min_degree_ordering` | Minimum degree heuristic |
| `ordering_tools` | Ordering utility functions |
| `reductions` | Data reduction rules (degree-1, simplicial) |

### `lib/spac/`

SPAC (Split-and-Connect) edge partitioning algorithm. Converts a node-partitioned graph into an edge-partitioned one using a bipartite graph construction.

### `lib/parallel_mh/`

Parallel memetic algorithm (used by kaffpaE):

| File | Description |
|---|---|
| `parallel_mh_async` | Asynchronous evolutionary loop |
| `population` | Manages population of partitions |
| `galinier_combine/` | Galinier combine operator (partition combination) |
| `exchange/` | MPI-based partition exchange between processes |

### `parallel/parallel_src/`

ParHIP — the distributed-memory parallel partitioner:

```
parallel/parallel_src/
├── interface/
│   ├── parhip_interface.h      # Public ParHIP C API
│   └── parhip_interface.cpp    # API implementation
├── app/
│   ├── parhip.cpp              # parhip executable
│   ├── toolbox.cpp             # toolbox executable
│   ├── graph2binary.cpp        # format converter
│   └── dspac.cpp               # distributed edge partitioning
└── lib/
    ├── distributed_partitioning/  # Distributed coarsening + refinement
    ├── parallel_label_compress/   # Parallel label compression
    └── tools/                     # Distributed graph I/O, metrics
```

---

## Data Flow: Serial Partitioning

```
User calls kaffpa()
    │
    ▼
kaHIP_interface.cpp
    │  Builds PartitionConfig from mode (FAST/ECO/STRONG)
    │  Reads graph data (CSR arrays) into graph_access
    ▼
graph_partitioner::perform_partitioning()
    │
    ├─── If k == 2: recursive_bisection()
    │
    └─── If k > 2: multilevel_recursion()
              │
              ▼
         wcycle_partitioner::perform_partitioning()
              │
              ├─1── coarsening::perform_coarsening()
              │       │  Loops until graph is small enough:
              │       │    1. edge_ratings: score each edge
              │       │    2. matching: find matching/clustering
              │       │    3. contraction: merge matched nodes
              │       ▼
              │     graph_hierarchy (stack of coarsened graphs)
              │
              ├─2── initial_partitioning::perform_initial_partitioning()
              │       │  At coarsest level (≈ 2k nodes):
              │       │    Multiple attempts of bipartition/recursion
              │       │    Keep best result
              │       ▼
              │     initial partition on coarsest graph
              │
              └─3── uncoarsening::perform_uncoarsening()
                      │  For each level in graph_hierarchy (bottom up):
                      │    1. Project partition to finer graph
                      │    2. refinement::perform_refinement()
                      │         → 2-way FM + max-flow (quotient graph)
                      │         → k-way FM (full graph)
                      │         → Label propagation
                      │         → Cycle improvements (if KaBaPE)
                      ▼
                    final partition (part[] array)
```

---

## Data Flow: Parallel Partitioning (ParHIP)

```
Each MPI rank calls ParHIPPartitionKWay()
    │
    ▼
Distributed graph stored as vtxdist/xadj/adjncy across ranks
    │
    ├─1── Distributed coarsening
    │       Label propagation with size constraints
    │       Parallel contraction
    │
    ├─2── Gather coarsest graph on rank 0
    │       Initial partitioning using full KaFFPa
    │       (or parallel memetic algorithm for high quality)
    │
    ├─3── Broadcast and project partition
    │
    └─4── Distributed refinement
            Label propagation with boundary-aware moves
```

---

## Configuration System

`PartitionConfig` (in `lib/partition/partition_config.h`) is a plain struct with 100+ fields. Preconfigurations (FAST, ECO, STRONG, …) fill this struct with tuned values. Users can override individual fields through the CLI or by modifying the struct programmatically.

```cpp
// Simplified flow inside kaHIP_interface.cpp:
PartitionConfig config;
switch(mode) {
    case FAST:   configuration.fast(config);   break;
    case ECO:    configuration.eco(config);    break;
    case STRONG: configuration.strong(config); break;
    // ...
}
config.k = nparts;
config.imbalance = imbalance * 100;  // convert from fraction to percent
```

---

## CMake Object Library Design

KaHIP uses CMake OBJECT libraries to avoid linking the same code multiple times into different executables. Each component is compiled once and its object files are linked directly:

```cmake
add_library(libkaffpa OBJECT ${LIBKAFFPA_SOURCE_FILES})
add_executable(kaffpa app/kaffpa.cpp $<TARGET_OBJECTS:libkaffpa>)
add_library(kahip SHARED interface/kaHIP_interface.cpp
                          $<TARGET_OBJECTS:libkaffpa>
                          $<TARGET_OBJECTS:libmapping>
                          $<TARGET_OBJECTS:libnodeordering>
                          $<TARGET_OBJECTS:libspac>)
```

This means `libkahip.so` is a true self-contained shared library including all of KaHIP's functionality.

---

## Thread Safety

- Serial KaFFPa functions (`kaffpa`, `kaffpa_balance`, etc.) are **not thread-safe** due to global state in the PRNG and configuration system.
- Each call to `kaffpa()` should be made from a single thread, or with external synchronization.
- OpenMP is used internally for parallel loops within a single partitioning call; this is safe as long as calls are serialized at the application level.
- ParHIP is inherently multi-process (MPI) and thread-safe within each rank.

---

## Index Type System

KaHIP supports both 32-bit and 64-bit edge indices via a compile-time typedef:

```c
// interface/kaHIP_interface.h
#ifdef KAHIP_64BIT
typedef int64_t kahip_idx;  // compile with -D64BITMODE=ON
#else
typedef int32_t kahip_idx;  // default
#endif
```

Node-related parameters (`n`, `vwgt`, `nparts`, `part`) always use `int` (32-bit).

Edge-related parameters (`xadj`, `adjncy`, `adjcwgt`, `edgecut`) use `kahip_idx`.

Use `kahip_sizeof_idx()` at runtime to query which mode was compiled in.
