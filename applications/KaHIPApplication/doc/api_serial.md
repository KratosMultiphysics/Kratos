# Serial C/C++ API Reference

Header: `interface/kaHIP_interface.h`
Library: `libkahip.so` / `libkahip_static.a`

---

## Index Type

```c
#ifdef KAHIP_64BIT
typedef int64_t kahip_idx;   // compile with -D64BITMODE=ON
#else
typedef int32_t kahip_idx;   // default (32-bit)
#endif

int kahip_sizeof_idx();      // returns 4 (32-bit) or 8 (64-bit) at runtime
```

Edge-related arrays (`xadj`, `adjncy`, `adjcwgt`) and scalar values (`edgecut`) use `kahip_idx`.
Node-related parameters (`n`, `vwgt`, `nparts`, `part`) always use `int`.

---

## Mode Constants

### Partitioning Modes

| Constant       | Value | Description                                  |
|----------------|-------|----------------------------------------------|
| `FAST`         | 0     | Speed-optimized; good for meshes             |
| `ECO`          | 1     | Balanced quality/speed; good default         |
| `STRONG`       | 2     | High quality; slower                         |
| `FASTSOCIAL`   | 3     | Speed-optimized for social/scale-free graphs |
| `ECOSOCIAL`    | 4     | Balanced for social graphs                   |
| `STRONGSOCIAL` | 5     | High quality for social graphs               |

### Mapping Modes

| Constant               | Value | Description                         |
|------------------------|-------|-------------------------------------|
| `MAPMODE_MULTISECTION` | 0     | Global multisection process mapping |
| `MAPMODE_BISECTION`    | 1     | Bisection-based process mapping     |

---

## Graph Representation

All functions use **CSR (Compressed Sparse Row)** format. For an undirected graph with `n` nodes and `m` edges:

```
n               : number of nodes (vertices)
xadj[0..n]      : xadj[i] is the start of node i's adjacency list in adjncy
                  xadj[n] = 2*m (total number of directed edge entries)
adjncy[0..2m-1] : neighbor node IDs (0-indexed)
vwgt[0..n-1]    : node weights; NULL means all weights = 1
adjcwgt[0..2m-1]: edge weights; NULL means all weights = 1
```

**Important**: Each undirected edge (u, v) appears twice: once in u's list and once in v's list. Edge weights must be consistent: `adjcwgt` at position of (u→v) must equal `adjcwgt` at position of (v→u).

Example — a path graph with 5 nodes (0—1—2—3—4):
```
n = 5
xadj    = [0, 1, 3, 5, 7, 8]
adjncy  = [1, 0, 2, 1, 3, 2, 4, 3]
vwgt    = NULL  (uniform)
adjcwgt = NULL  (uniform)
```

---

## Functions

### `kaffpa` — Main Partitioning

```c
void kaffpa(
    int*        n,              // [in]  number of nodes
    int*        vwgt,           // [in]  node weights array (NULL = uniform)
    kahip_idx*  xadj,           // [in]  CSR row offsets (size n+1)
    kahip_idx*  adjcwgt,        // [in]  edge weights (NULL = uniform)
    kahip_idx*  adjncy,         // [in]  CSR column indices
    int*        nparts,         // [in]  number of blocks k
    double*     imbalance,      // [in]  allowed imbalance ε (e.g. 0.03 = 3%)
    bool        suppress_output,// [in]  suppress stdout if true
    int         seed,           // [in]  random seed (0 = default)
    int         mode,           // [in]  FAST / ECO / STRONG / *SOCIAL
    kahip_idx*  edgecut,        // [out] total edge cut weight
    int*        part            // [out] partition array (size n), pre-allocated
);
```

**Description**: Partitions the graph into `nparts` blocks. The partition minimizes the weighted edge cut subject to the balance constraint: each block has weight at most `(1 + imbalance) * total_weight / nparts`.

**Notes**:
- `part` must be pre-allocated by the caller with `n` integers.
- `edgecut` is the sum of weights of edges crossing block boundaries.
- Setting `seed = 0` uses a fixed default seed for reproducibility.
- Increasing the seed value changes the random choices; running multiple seeds and taking the best result improves quality.

**Example**:
```cpp
#include "kaHIP_interface.h"

int n = 5;
kahip_idx xadj[]   = {0, 2, 5, 7, 9, 12};
kahip_idx adjncy[] = {1,4, 0,2,4, 1,3, 2,4, 0,1,3};

int nparts         = 2;
double imbalance   = 0.03;
int* part          = new int[n];
kahip_idx edgecut  = 0;

kaffpa(&n, nullptr, xadj, nullptr, adjncy,
       &nparts, &imbalance, false, 0, ECO,
       &edgecut, part);

// part[i] ∈ {0, 1, …, nparts-1}
```

---

### `kaffpa_balance` — Partitioning with Perfect Balance Option

```c
void kaffpa_balance(
    int*        n,
    int*        vwgt,
    kahip_idx*  xadj,
    kahip_idx*  adjcwgt,
    kahip_idx*  adjncy,
    int*        nparts,
    double*     imbalance,
    bool        perfectly_balance,  // [in] enforce exact balance if true
    bool        suppress_output,
    int         seed,
    int         mode,
    kahip_idx*  edgecut,
    int*        part
);
```

**Description**: Like `kaffpa`, but with an optional `perfectly_balance` flag. When `perfectly_balance = true`, the algorithm uses the KaBaPE approach to enforce `w(Vᵢ) = W/k` exactly for all blocks (or as close as integer weights allow). This may increase edge cut.

---

### `kaffpa_balance_NE` — Node and Edge Balance

```c
void kaffpa_balance_NE(
    int*        n,
    int*        vwgt,
    kahip_idx*  xadj,
    kahip_idx*  adjcwgt,
    kahip_idx*  adjncy,
    int*        nparts,
    double*     imbalance,
    bool        suppress_output,
    int         seed,
    int         mode,
    kahip_idx*  edgecut,
    int*        part
);
```

**Description**: Enforces balance constraints on both **node weights** and **edge weights** simultaneously. Useful when both node computational load and communication volume per block must be balanced.

---

### `node_separator` — Vertex Separator

```c
void node_separator(
    int*        n,
    int*        vwgt,
    kahip_idx*  xadj,
    kahip_idx*  adjcwgt,
    kahip_idx*  adjncy,
    int*        nparts,
    double*     imbalance,
    bool        suppress_output,
    int         seed,
    int         mode,
    int*        num_separator_vertices,  // [out] size of separator
    int**       separator                // [out] array of separator node IDs
);
```

**Description**: Computes a node separator S such that removing S from the graph leaves no edges between the remaining parts. The result partitions V into (A, B, S) where A and B are disconnected.

For k > 2, the separator connects all k parts. The algorithm first computes a k-way partition, then uses flow-based local search to find an optimal separator.

**Memory**: `*separator` is allocated by KaHIP and must be freed by the caller with `free(*separator)`.

**Example**:
```cpp
int num_sep;
int* separator;
node_separator(&n, nullptr, xadj, nullptr, adjncy,
               &nparts, &imbalance, false, 0, STRONG,
               &num_sep, &separator);
// Use separator[0..num_sep-1]
free(separator);
```

---

### `process_mapping` — QAP-Based Process Mapping

```c
void process_mapping(
    int*        n,
    int*        vwgt,
    kahip_idx*  xadj,
    kahip_idx*  adjcwgt,
    kahip_idx*  adjncy,
    int*        hierarchy_parameter,   // [in] processor hierarchy sizes
    int*        distance_parameter,    // [in] communication distances
    int         hierarchy_depth,       // [in] depth of hierarchy
    int         mode_partitioning,     // [in] FAST / ECO / STRONG
    int         mode_mapping,          // [in] MAPMODE_MULTISECTION / MAPMODE_BISECTION
    double*     imbalance,
    bool        suppress_output,
    int         seed,
    kahip_idx*  edgecut,               // [out] edge cut
    int*        qap,                   // [out] QAP objective value
    int*        part                   // [out] partition/mapping array
);
```

**Description**: Maps the graph nodes to processors in a hierarchically organized system. The `hierarchy_parameter` array describes the system topology (e.g., `[4, 8, 8]` = 4 nodes × 8 sockets × 8 cores = 256 PEs). The `distance_parameter` array gives the communication cost at each level (e.g., `[1, 10, 100]`).

**Example** (4×8×8 = 256-core system):
```cpp
int hierarchy[] = {4, 8, 8};
int distances[] = {1, 10, 100};
int qap_val;
process_mapping(&n, vwgt, xadj, adjcwgt, adjncy,
                hierarchy, distances, 3,
                ECO, MAPMODE_MULTISECTION,
                &imbalance, false, 0,
                &edgecut, &qap_val, part);
```

---

### `reduced_nd` — Node Ordering (Nested Dissection)

```c
void reduced_nd(
    int*        n,
    kahip_idx*  xadj,
    kahip_idx*  adjncy,
    bool        suppress_output,
    int         seed,
    int         mode,           // use STRONG for best fill-in reduction
    int*        ordering        // [out] permutation array (size n)
);
```

**Description**: Computes a fill-in reducing node ordering using data reduction rules followed by multilevel nested dissection. The `ordering` array is a permutation of `{0, 1, …, n-1}`. If the graph is interpreted as a sparse matrix A, then `A[ordering, ordering]` has reduced fill-in during factorization.

**Note**: The graph must be unweighted (edge weights ignored).

If compiled with METIS (`-DUSEMETIS`):
```c
void reduced_nd_fast(int* n, kahip_idx* xadj, kahip_idx* adjncy,
                     bool suppress_output, int seed, int* ordering);
```
This applies reduction rules before calling METIS nested dissection — typically faster and better quality.

---

### `edge_partitioning` — Edge-Centric Partitioning

```c
void edge_partitioning(
    int*        n,
    int*        vwgt,
    kahip_idx*  xadj,
    kahip_idx*  adjcwgt,
    kahip_idx*  adjncy,
    int*        nparts,
    double*     imbalance,
    bool        suppress_output,
    int         seed,
    int         mode,
    int*        vertexcut,              // [out] vertex replication count
    int*        part,                   // [out] edge partition labels (size = number of edges)
    kahip_idx   infinity_edge_weight    // [in]  weight for infinity edges (default 1000)
);
```

**Description**: Partitions the **edges** of the graph (not nodes) into k equally-sized blocks, minimizing vertex replications. This is useful for "think-like-an-edge" distributed computation frameworks.

The output `part` array has one entry per **directed edge** in the adjacency list (size = `xadj[n]`).

---

### `kahip_sizeof_idx` — Runtime Index Size Query

```c
int kahip_sizeof_idx();
```

Returns `4` if compiled in 32-bit mode (default), `8` if compiled with `-D64BITMODE=ON`.

Use this to safely query the active mode when using the library from C or dynamically linked contexts.

---

## Linking

### CMake

```cmake
find_package(PkgConfig REQUIRED)
pkg_check_modules(KAHIP REQUIRED kahip)

target_include_directories(myapp PRIVATE ${KAHIP_INCLUDE_DIRS})
target_link_libraries(myapp PRIVATE ${KAHIP_LIBRARIES})
```

### Manual

```bash
# Compile
g++ myapp.cpp -I/path/to/kahip/include -L/path/to/kahip/lib -lkahip -fopenmp -o myapp

# Or with static library
g++ myapp.cpp -I/path/to/kahip/include /path/to/kahip/lib/libkahip.a -fopenmp -o myapp
```

### pkg-config

```bash
g++ myapp.cpp $(pkg-config --cflags --libs kahip) -o myapp
```

---

## Complete C++ Example

```cpp
#include <iostream>
#include <vector>
#include "kaHIP_interface.h"

int main() {
    // Triangle graph: 0-1, 1-2, 2-0
    int n = 3;
    std::vector<kahip_idx> xadj   = {0, 2, 4, 6};
    std::vector<kahip_idx> adjncy = {1, 2,  0, 2,  0, 1};
    std::vector<int>       vwgt   = {1, 1, 1};
    std::vector<kahip_idx> adjcwgt = {1, 1, 1, 1, 1, 1};

    int nparts       = 2;
    double imbalance = 0.03;
    std::vector<int> part(n);
    kahip_idx edgecut = 0;

    kaffpa(&n,
           vwgt.data(),
           xadj.data(),
           adjcwgt.data(),
           adjncy.data(),
           &nparts,
           &imbalance,
           /*suppress_output=*/true,
           /*seed=*/42,
           ECO,
           &edgecut,
           part.data());

    std::cout << "Edge cut: " << edgecut << "\n";
    for (int i = 0; i < n; ++i)
        std::cout << "Node " << i << " → block " << part[i] << "\n";

    return 0;
}
```

---

## Complete C Example

```c
#include <stdio.h>
#include <stdlib.h>
#include "kaHIP_interface.h"

int main(void) {
    int n = 5;
    kahip_idx xadj[]   = {0, 2, 5, 7, 9, 12};
    kahip_idx adjncy[] = {1,4, 0,2,4, 1,3, 2,4, 0,1,3};

    int nparts = 2;
    double imbalance = 0.03;
    int* part = (int*)malloc(n * sizeof(int));
    kahip_idx edgecut = 0;

    kaffpa(&n, NULL, xadj, NULL, adjncy,
           &nparts, &imbalance, 1, 0, ECO,
           &edgecut, part);

    printf("Edge cut: %d\n", (int)edgecut);
    for (int i = 0; i < n; i++)
        printf("Node %d -> block %d\n", i, part[i]);

    free(part);
    return 0;
}
```
