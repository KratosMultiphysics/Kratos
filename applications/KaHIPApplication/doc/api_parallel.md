# Parallel API Reference

KaHIP provides two distinct parallel partitioning approaches:

1. **ParHIP** (`parhip_interface.h`) — distributed-memory MPI-based parallel partitioner designed for large-scale complex networks
2. **kaffpaE** — MPI-parallel evolutionary (memetic) algorithm that wraps the serial KaFFPa for high-quality solutions

---

## ParHIP API

Header: `parallel/parallel_src/interface/parhip_interface.h`
Library: `libparhip_interface.a` / `libparhip_interface.so`

### Index Type

```c
typedef unsigned long long idxtype;   // always 64-bit in ParHIP
```

Note: ParHIP uses `unsigned long long` (`uint64_t`) for all indices, unlike the serial API which uses `kahip_idx`.

### Mode Constants

| Constant | Value | Description |
|---|---|---|
| `ULTRAFASTMESH` | 0 | Fastest, for mesh graphs |
| `FASTMESH` | 1 | Fast, for mesh graphs |
| `ECOMESH` | 2 | Balanced quality for meshes |
| `ULTRAFASTSOCIAL` | 3 | Fastest, for social/complex networks |
| `FASTSOCIAL` | 4 | Fast, for social/complex networks |
| `ECOSOCIAL` | 5 | Balanced quality for social networks |

### `ParHIPPartitionKWay`

```c
void ParHIPPartitionKWay(
    idxtype*    vtxdist,        // [in]  distribution of vertices across ranks
    idxtype*    xadj,           // [in]  local CSR row offsets (size local_n+1)
    idxtype*    adjncy,         // [in]  local CSR column indices
    idxtype*    vwgt,           // [in]  local node weights (NULL = uniform)
    idxtype*    adjwgt,         // [in]  local edge weights (NULL = uniform)
    int*        nparts,         // [in]  total number of blocks k
    double*     imbalance,      // [in]  allowed imbalance ε (e.g. 0.03)
    bool        suppress_output,// [in]  suppress stdout
    int         seed,           // [in]  random seed
    int         mode,           // [in]  ULTRAFASTMESH, FASTMESH, ECOMESH, ULTRAFASTSOCIAL, FASTSOCIAL, ECOSOCIAL
    int*        edgecut,        // [out] global edge cut
    idxtype*    part,           // [out] local partition array (size local_n)
    MPI_Comm*   comm            // [in]  MPI communicator
);
```

### Distributed Graph Format

ParHIP uses a distributed CSR format. The graph is split across all MPI ranks. The `vtxdist` array describes this split:

```
vtxdist[p]   = first global node ID owned by rank p
vtxdist[p+1] = first global node ID owned by rank p+1

Rank p owns nodes [vtxdist[p], vtxdist[p+1])
local_n[p] = vtxdist[p+1] - vtxdist[p]
```

`xadj` and `adjncy` use **local** indices for the owned nodes, but **global** indices for neighbors (i.e., `adjncy` entries are global node IDs, 0-indexed).

**Example for 3 ranks, 9 nodes (3 per rank):**
```
vtxdist = [0, 3, 6, 9]

Rank 0: owns nodes 0,1,2; xadj[0..3], adjncy has global IDs
Rank 1: owns nodes 3,4,5; xadj[0..3], adjncy has global IDs
Rank 2: owns nodes 6,7,8; xadj[0..3], adjncy has global IDs
```

### Minimal ParHIP Example

```cpp
#include <mpi.h>
#include "parhip_interface.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Total graph: 8 nodes arranged in a ring
    idxtype total_n = 8;

    // Distribute evenly: vtxdist[p] = p * (total_n / size)
    std::vector<idxtype> vtxdist(size + 1);
    for (int p = 0; p <= size; p++)
        vtxdist[p] = (idxtype)p * total_n / size;

    idxtype local_n = vtxdist[rank + 1] - vtxdist[rank];

    // Build local portion of ring adjacency (each node connects to prev/next)
    std::vector<idxtype> xadj(local_n + 1);
    std::vector<idxtype> adjncy;
    xadj[0] = 0;
    for (idxtype i = 0; i < local_n; i++) {
        idxtype global_i = vtxdist[rank] + i;
        adjncy.push_back((global_i + total_n - 1) % total_n);  // prev
        adjncy.push_back((global_i + 1) % total_n);             // next
        xadj[i + 1] = adjncy.size();
    }

    int nparts = 2;
    double imbalance = 0.03;
    int edgecut = 0;
    std::vector<idxtype> part(local_n);

    ParHIPPartitionKWay(vtxdist.data(), xadj.data(), adjncy.data(),
                        nullptr, nullptr,
                        &nparts, &imbalance,
                        /*suppress_output=*/false, /*seed=*/0,
                        FASTMESH, &edgecut, part.data(), &comm);

    if (rank == 0)
        printf("Edge cut: %d\n", edgecut);

    MPI_Finalize();
    return 0;
}
```

Compile and run:
```bash
mpicxx -o parhip_test parhip_test.cpp \
    -I/path/to/kahip/include \
    -L/path/to/kahip/lib -lparhip_interface

mpirun -np 4 ./parhip_test
```

### ParHIP-specific Notes

- **Input format**: ParHIP can read METIS `.graph` text files or binary `.bgf` files. Binary files are much faster for large graphs.
- **Binary conversion**: Use `graph2binary` to convert METIS format to BGF before distributing.
- **Determinism**: Use `-DDETERMINISTIC_PARHIP=ON` at build time for reproducible results (reduces quality slightly).
- **Output**: The full partition is assembled on rank 0; other ranks receive their local portion.

---

## kaffpaE — Parallel Evolutionary Algorithm

kaffpaE is an MPI-parallel application (not a library function) that runs a parallel memetic (evolutionary) algorithm. It uses KaFFPa as the mutation/combine operator and exchanges partitions between MPI ranks to improve quality over time.

### Running kaffpaE

```bash
# Basic usage: 24 MPI ranks, time limit 3600s
mpirun -n 24 ./kaffpaE graph.metis --k 4 --time_limit=3600

# High quality with all enhancements
mpirun -n 24 ./kaffpaE graph.metis \
    --k 4 \
    --time_limit=3600 \
    --preconfiguration=strong \
    --mh_enable_tabu_search \
    --mh_enable_kabapE

# Social network mode
mpirun -n 24 ./kaffpaE graph.metis \
    --k 32 \
    --time_limit=7200 \
    --preconfiguration=ssocial

# With connected blocks (experimental, requires connected input)
mpirun -n 24 ./kaffpaE graph.metis \
    --k 8 \
    --time_limit=1800 \
    --preconfiguration=strong \
    --connected_blocks
```

### Key kaffpaE Options

| Option | Default | Description |
|---|---|---|
| `--k` | required | Number of blocks |
| `--time_limit` | 0 | Wall-clock time limit in seconds |
| `--preconfiguration` | `strong` | `fast`, `eco`, `strong`, `fsocial`, `esocial`, `ssocial` |
| `--seed` | 0 | Random seed |
| `--mh_enable_tabu_search` | off | Enable tabu search improvement |
| `--mh_enable_kabapE` | off | Enable KaBaPE combine operator |
| `--mh_enable_gal_combine` | off | Enable Galinier combine operator |
| `--imbalance` | 3 | Allowed imbalance in percent |
| `--suppress_output` | off | Suppress verbose output |
| `--output_filename` | — | Write partition to file |
| `--connected_blocks` | off | Enforce connected partition blocks (strong only) |

### Algorithm

kaffpaE runs an asynchronous parallel memetic algorithm:

1. Each MPI rank starts with an independent partition (from KaFFPa)
2. Periodically, ranks exchange their best partitions
3. The combine operator (Galinier or KaBaPE) merges two parent partitions
4. KaFFPa is used as a mutation operator
5. The best partition found across all ranks and all time is returned

This is significantly higher quality than serial KaFFPa, especially for large `k` and complex graphs.

---

## ParHIP CLI Programs

### `parhip`

Distributed parallel graph partitioner.

```bash
mpirun -n <P> ./parhip graph.graph --k <k> --preconfiguration=<mode>

# Modes: ultrafastmesh, fastmesh, ecomesh, ultrafastsocial, fastsocial, ecosocial
# Options:
#   --k             number of blocks (required)
#   --seed          random seed
#   --imbalance     allowed imbalance in percent
#   --output_filename FILE   write partition to FILE
#   --suppress_output        suppress output
```

```bash
# Fast mesh partitioning
mpirun -n 8 ./parhip rgg_n_2_15_s0.graph --k 16 --preconfiguration=fastmesh

# Social network (binary input for speed)
./graph2binary social.graph social.bgf
mpirun -n 32 ./parhip social.bgf --k 64 --preconfiguration=ecosocial
```

### `graph2binary`

Convert METIS format graph to binary BGF format for faster ParHIP I/O.

```bash
./graph2binary input.graph output.bgf
```

### `graph2binary_external`

Convert a graph with an externally provided partition.

```bash
./graph2binary_external input.graph input.partition output.bgf
```

### `toolbox`

Evaluate partition quality and convert between formats.

```bash
./toolbox graph.graph partition.txt --evaluate
```

### `dspac` (Distributed SPAC)

Distributed edge partitioning.

```bash
mpirun -n 8 ./distributed_edge_partitioning graph.bgf --k 16 --preconfiguration=fastsocial
```

---

## Choosing Between ParHIP and kaffpaE

| Aspect | ParHIP | kaffpaE |
|---|---|---|
| **Input** | Distributed across all ranks | Each rank reads full file |
| **Memory** | Distributed (scales to huge graphs) | Full graph on each rank |
| **Graph type** | Social/complex networks (best), meshes (ok) | Any graph type |
| **Quality** | Good (label propagation based) | Excellent (evolutionary) |
| **Speed** | Fastest for large graphs | Slower but better quality |
| **Scalability** | Excellent (distributed memory) | Limited by single-node memory |
| **When to use** | Graphs too large for one node | High quality, fits in RAM |

**Rule of thumb**: Use ParHIP for graphs with billions of edges (distributed memory required). Use kaffpaE when the graph fits in RAM but you need higher quality than serial KaFFPa.

---

## Linking ParHIP

### CMake

```cmake
find_package(PkgConfig REQUIRED)
pkg_check_modules(PARHIP REQUIRED parhip)
find_package(MPI REQUIRED)

target_include_directories(myapp PRIVATE ${PARHIP_INCLUDE_DIRS})
target_link_libraries(myapp PRIVATE ${PARHIP_LIBRARIES} MPI::MPI_CXX)
```

### Manual

```bash
mpicxx myapp.cpp \
    -I/path/to/kahip/include \
    -L/path/to/kahip/lib \
    -lparhip_interface \
    -o myapp

mpirun -np 4 ./myapp
```
