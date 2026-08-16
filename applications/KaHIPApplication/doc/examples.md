# Usage Examples

## CLI Examples

### Validate a graph file

```bash
./graphchecker examples/rgg_n_2_15_s0.graph
```

### Serial partitioning — kaffpa

```bash
# Fast partitioning into 4 blocks
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --preconfiguration=fast

# Eco quality (good default)
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --preconfiguration=eco

# Strong (highest serial quality)
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --preconfiguration=strong

# With imbalance tolerance of 5%
./kaffpa examples/rgg_n_2_15_s0.graph --k 8 --preconfiguration=eco --imbalance=5

# Save partition to file
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --output_filename=partition.txt

# Social network mode
./kaffpa social_graph.graph --k 16 --preconfiguration=esocial

# Start from existing partition
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --input_partition=partition.txt --preconfiguration=strong

# With specific seed for reproducibility
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --seed=12345 --preconfiguration=eco

# Connected blocks (experimental, input must be connected)
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --preconfiguration=strong --connected_blocks

# Fast I/O for large files
./kaffpa large_graph.graph --k 64 --preconfiguration=fast --mmap_io

# Process mapping (4 nodes × 8 sockets × 8 cores = 256 PEs)
./kaffpa examples/rgg_n_2_15_s0.graph --k 256 \
    --preconfiguration=eco \
    --enable_mapping \
    --hierarchy_parameter_string=4:8:8 \
    --distance_parameter_string=1:10:100
```

### Parallel evolutionary partitioning — kaffpaE

```bash
# High quality with 24 MPI processes, 1 hour time limit
mpirun -n 24 ./kaffpaE examples/rgg_n_2_15_s0.graph \
    --k 4 \
    --time_limit=3600 \
    --preconfiguration=strong

# Maximum quality with all enhancement options
mpirun -n 24 ./kaffpaE examples/rgg_n_2_15_s0.graph \
    --k 4 \
    --time_limit=3600 \
    --preconfiguration=strong \
    --mh_enable_tabu_search \
    --mh_enable_kabapE

# Social network with eco social preconfiguration
mpirun -n 16 ./kaffpaE social_graph.graph \
    --k 32 \
    --time_limit=7200 \
    --preconfiguration=ssocial

# Connected blocks (experimental)
mpirun -n 8 ./kaffpaE examples/rgg_n_2_15_s0.graph \
    --k 8 \
    --time_limit=600 \
    --preconfiguration=strong \
    --connected_blocks \
    --output_filename=partition_connected.txt
```

### Distributed parallel partitioning — ParHIP

```bash
# Convert text format to binary (much faster I/O for ParHIP)
./graph2binary examples/rgg_n_2_15_s0.graph examples/rgg_n_2_15_s0.bgf

# Fast mesh partitioning
mpirun -n 8 ./parhip examples/rgg_n_2_15_s0.bgf \
    --k 32 \
    --preconfiguration=fastmesh

# Eco mesh quality
mpirun -n 16 ./parhip examples/rgg_n_2_15_s0.graph \
    --k 32 \
    --preconfiguration=ecomesh

# Social network (ultra-fast mode)
mpirun -n 32 ./parhip large_social.bgf \
    --k 128 \
    --preconfiguration=ultrafastsocial

# Save partition
mpirun -n 8 ./parhip examples/rgg_n_2_15_s0.bgf \
    --k 16 \
    --preconfiguration=fastmesh \
    --output_filename=partition_parallel.txt
```

### Node separator

```bash
# 2-way node separator
./node_separator examples/rgg_n_2_15_s0.graph

# k-way separator: first partition, then extract separator
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --output_filename=partition.txt
./partition_to_vertex_separator examples/rgg_n_2_15_s0.graph --k 4 --input_partition=partition.txt
```

### Node ordering (sparse matrix fill-in reduction)

```bash
# Standard nested dissection ordering
./node_ordering examples/rgg_n_2_15_s0.graph

# Fast ordering using reductions + METIS (if METIS installed)
./fast_node_ordering examples/rgg_n_2_15_s0.graph
```

### Edge partitioning

```bash
# Sequential edge partitioning
./edge_partitioning examples/rgg_n_2_15_s0.graph --k 4 --preconfiguration=fast

# Distributed edge partitioning (ParHIP required)
mpirun -n 4 ./distributed_edge_partitioning examples/rgg_n_2_15_s0.bgf \
    --k 4 \
    --preconfiguration=fastsocial
```

### Partition evaluation

```bash
# Evaluate quality of an existing partition
./evaluator examples/rgg_n_2_15_s0.graph --k 4 --input_partition=partition.txt

# Label propagation clustering
./label_propagation examples/rgg_n_2_15_s0.graph
```

---

## C++ Library Examples

### Minimal C++ Example

```cpp
#include <iostream>
#include "kaHIP_interface.h"

int main() {
    // Triangle: nodes 0,1,2; edges 0-1, 1-2, 0-2
    int n = 3;
    kahip_idx xadj[]   = {0, 2, 4, 6};
    kahip_idx adjncy[] = {1, 2, 0, 2, 0, 1};

    int nparts = 2;
    double imbalance = 0.03;
    int part[3];
    kahip_idx edgecut = 0;

    kaffpa(&n, nullptr, xadj, nullptr, adjncy,
           &nparts, &imbalance, true, 0, ECO,
           &edgecut, part);

    std::cout << "Edge cut: " << edgecut << "\n";
    return 0;
}
```

Compile:
```bash
g++ example.cpp -I/path/to/kahip/include -L/path/to/kahip/lib -lkahip -fopenmp -o example
```

### Loading METIS Graph and Partitioning

```cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "kaHIP_interface.h"

// Returns false on failure
bool load_metis_graph(const std::string& filename,
                      int& n,
                      std::vector<kahip_idx>& xadj,
                      std::vector<kahip_idx>& adjncy,
                      std::vector<int>& vwgt,
                      std::vector<kahip_idx>& adjcwgt) {
    std::ifstream f(filename);
    if (!f) return false;

    std::string line;
    while (std::getline(f, line))
        if (line[0] != '%') break;

    std::istringstream hdr(line);
    int m, fmt = 0;
    hdr >> n >> m >> fmt;

    bool has_vwgt = (fmt / 10) % 10 == 1;
    bool has_ewgt = fmt % 10 == 1;

    xadj.resize(n + 1);
    vwgt.resize(n, 1);
    adjncy.reserve(2 * m);
    if (has_ewgt) adjcwgt.reserve(2 * m);

    xadj[0] = 0;
    for (int i = 0; i < n; ++i) {
        while (std::getline(f, line) && line[0] == '%') {}
        std::istringstream row(line);
        if (has_vwgt) row >> vwgt[i];
        int neighbor;
        while (row >> neighbor) {
            adjncy.push_back(neighbor - 1);  // 1-indexed → 0-indexed
            if (has_ewgt) {
                int ew; row >> ew;
                adjcwgt.push_back(ew);
            } else {
                adjcwgt.push_back(1);
            }
        }
        xadj[i + 1] = adjncy.size();
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <graph.graph> <k>\n";
        return 1;
    }

    int n;
    std::vector<kahip_idx> xadj, adjncy, adjcwgt;
    std::vector<int> vwgt;

    if (!load_metis_graph(argv[1], n, xadj, adjncy, vwgt, adjcwgt)) {
        std::cerr << "Failed to load graph\n";
        return 1;
    }

    int nparts = std::stoi(argv[2]);
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
           false, 0, ECO,
           &edgecut,
           part.data());

    std::cout << "Nodes: " << n << ", Parts: " << nparts << "\n";
    std::cout << "Edge cut: " << edgecut << "\n";

    // Print block sizes
    std::vector<int> block_sizes(nparts, 0);
    for (int i = 0; i < n; ++i) block_sizes[part[i]]++;
    for (int k = 0; k < nparts; ++k)
        std::cout << "Block " << k << ": " << block_sizes[k] << " nodes\n";

    return 0;
}
```

### Trying Multiple Seeds

```cpp
#include "kaHIP_interface.h"
#include <vector>
#include <limits>

std::pair<kahip_idx, std::vector<int>> partition_best_of_n(
    int n, std::vector<kahip_idx>& xadj, std::vector<kahip_idx>& adjncy,
    int nparts, double imbalance, int n_tries = 10, int base_mode = ECO)
{
    kahip_idx best_cut = std::numeric_limits<kahip_idx>::max();
    std::vector<int> best_part(n), part(n);

    for (int seed = 0; seed < n_tries; ++seed) {
        kahip_idx cut = 0;
        kaffpa(&n, nullptr, xadj.data(), nullptr, adjncy.data(),
               &nparts, &imbalance, true, seed, base_mode, &cut, part.data());
        if (cut < best_cut) {
            best_cut = cut;
            best_part = part;
        }
    }
    return {best_cut, best_part};
}
```

---

## Python Examples

### Quick Partition

```python
import kahip

xadj   = [0, 2, 5, 7, 9, 12]
adjncy = [1, 4, 0, 2, 4, 1, 3, 2, 4, 0, 1, 3]
vwgt   = [1, 1, 1, 1, 1]
adjcwgt = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

edgecut, blocks = kahip.kaffpa(
    vwgt, xadj, adjcwgt, adjncy,
    nblocks=2, imbalance=0.03, suppress_output=0, seed=0, mode=kahip.ECO
)
print(f"Edge cut: {edgecut}, Blocks: {blocks}")
```

### Graph Builder Class

```python
import kahip

g = kahip.kahip_graph()
g.set_num_nodes(5)
g.add_undirected_edge(0, 1, 2)
g.add_undirected_edge(1, 2, 1)
g.add_undirected_edge(2, 3, 3)
g.add_undirected_edge(3, 4, 1)
g.add_undirected_edge(0, 4, 5)
g.set_weight(0, 2)

vwgt, xadj, adjcwgt, adjncy = g.get_csr_arrays()
edgecut, blocks = kahip.kaffpa(
    vwgt, xadj, adjcwgt, adjncy, 2, 0.03, 1, 0, kahip.STRONG
)
print(f"Edgecut: {edgecut}, Blocks: {blocks}")
```

---

## ParHIP C++ Library Example

```cpp
#include <mpi.h>
#include <iostream>
#include <vector>
#include "parhip_interface.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // 8-node ring graph, distributed across all ranks
    const idxtype total_n = 8;
    std::vector<idxtype> vtxdist(size + 1);
    for (int p = 0; p <= size; p++)
        vtxdist[p] = (idxtype)p * total_n / size;

    idxtype local_n = vtxdist[rank + 1] - vtxdist[rank];
    std::vector<idxtype> xadj(local_n + 1), adjncy;
    xadj[0] = 0;
    for (idxtype i = 0; i < local_n; i++) {
        idxtype g = vtxdist[rank] + i;
        adjncy.push_back((g + total_n - 1) % total_n);
        adjncy.push_back((g + 1) % total_n);
        xadj[i + 1] = adjncy.size();
    }

    int nparts = 2;
    double imbalance = 0.03;
    int edgecut = 0;
    std::vector<idxtype> part(local_n);

    ParHIPPartitionKWay(vtxdist.data(), xadj.data(), adjncy.data(),
                        nullptr, nullptr,
                        &nparts, &imbalance, false, 0,
                        FASTMESH, &edgecut, part.data(), &comm);

    if (rank == 0) printf("Edge cut: %d\n", edgecut);
    MPI_Finalize();
    return 0;
}
```

Compile:
```bash
mpicxx example_parhip.cpp \
    -I/path/to/kahip/include \
    -L/path/to/kahip/lib \
    -lparhip_interface \
    -o example_parhip

mpirun -np 4 ./example_parhip
```
