# Python Bindings

KaHIP provides Python bindings via pybind11, available both via PyPI and as a compiled module.

---

## Installation

### Via PyPI (recommended)

```bash
pip install kahip
```

This installs pre-built wheels for Python 3.9–3.14 on Linux, macOS, and Windows.

### From source (latest commit)

```bash
pip install pybind11

# From the KaHIP repository:
cd build/
./configure.sh --python

# Or manually:
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILDPYTHONMODULE=ON \
    -Dpybind11_DIR=$(python3 -m pip show pybind11 | grep Location | cut -d' ' -f2)/pybind11/share/cmake/pybind11
make -j$(nproc)

# The module is in build/kahip.cpython-<ver>-<platform>.so
# Add the build directory to PYTHONPATH to use it
export PYTHONPATH=$(pwd):$PYTHONPATH
```

### Homebrew (macOS/Linux)

```bash
brew install KaHIP/kahip/kahip
# Python bindings are included
```

---

## API Reference

### Module: `kahip`

```python
import kahip
```

### Mode Constants

```python
kahip.FAST        = 0
kahip.ECO         = 1
kahip.STRONG      = 2
kahip.FASTSOCIAL  = 3
kahip.ECOSOCIAL   = 4
kahip.STRONGSOCIAL = 5
```

### `kahip.kaffpa`

```python
edgecut, blocks = kahip.kaffpa(
    vwgt,           # list[int] or None — node weights (None = uniform)
    xadj,           # list[int] — CSR row offsets (length n+1)
    adjcwgt,        # list[int] or None — edge weights (None = uniform)
    adjncy,         # list[int] — CSR column indices
    nblocks,        # int — number of blocks k
    imbalance,      # float — allowed imbalance (e.g. 0.03 = 3%)
    suppress_output,# int — 0 = verbose, 1 = silent
    seed,           # int — random seed
    mode            # int — FAST, ECO, STRONG, etc.
)
# Returns:
#   edgecut: int — total weight of cut edges
#   blocks:  list[int] — block assignment for each node (0-indexed)
```

**Note**: The Python API argument order differs slightly from the C API.

### `kahip.kahip_graph`

Helper class for building graphs programmatically:

```python
g = kahip.kahip_graph()
g.set_num_nodes(n)              # set number of nodes
g.add_undirected_edge(u, v, w)  # add undirected edge (u,v) with weight w
g.set_weight(v, w)              # set node weight of node v to w
vwgt, xadj, adjcwgt, adjncy = g.get_csr_arrays()
```

---

## Examples

### Basic Partitioning

```python
import kahip

# Manual CSR representation (5-node path graph)
xadj   = [0, 2, 5, 7, 9, 12]
adjncy = [1, 4,  0, 2, 4,  1, 3,  2, 4,  0, 1, 3]
vwgt   = [1, 1, 1, 1, 1]
adjcwgt = [1, 1,  1, 1, 1,  1, 1,  1, 1,  1, 1, 1]

edgecut, blocks = kahip.kaffpa(
    vwgt, xadj, adjcwgt, adjncy,
    nblocks=2,
    imbalance=0.03,
    suppress_output=0,
    seed=42,
    mode=kahip.ECO
)

print(f"Edge cut: {edgecut}")
print(f"Node assignments: {blocks}")
# blocks[i] ∈ {0, 1, …, nblocks-1}
```

### Using `kahip_graph` Class

```python
import kahip

# Build a graph using the helper class
g = kahip.kahip_graph()
g.set_num_nodes(6)

# Add edges (undirected, with weights)
g.add_undirected_edge(0, 1, 1)
g.add_undirected_edge(1, 2, 1)
g.add_undirected_edge(2, 3, 1)
g.add_undirected_edge(3, 4, 1)
g.add_undirected_edge(4, 5, 1)
g.add_undirected_edge(0, 3, 10)  # cross-partition edge (heavy weight)

# Set node weights
for i in range(6):
    g.set_weight(i, 1)

# Get CSR and partition
vwgt, xadj, adjcwgt, adjncy = g.get_csr_arrays()

edgecut, blocks = kahip.kaffpa(
    vwgt, xadj, adjcwgt, adjncy,
    nblocks=2,
    imbalance=0.03,
    suppress_output=1,
    seed=0,
    mode=kahip.STRONG
)

print(f"Edgecut: {edgecut}")
print(f"Partition: {blocks}")
```

### High-Quality Partitioning with Multiple Seeds

```python
import kahip

def partition_best_of_n(xadj, adjncy, vwgt, adjcwgt, k, imbalance, n_trials=10, mode=kahip.ECO):
    """Try n_trials different seeds, return best partition."""
    best_cut = float('inf')
    best_blocks = None
    for seed in range(n_trials):
        cut, blocks = kahip.kaffpa(
            vwgt, xadj, adjcwgt, adjncy, k, imbalance,
            suppress_output=1, seed=seed, mode=mode
        )
        if cut < best_cut:
            best_cut = cut
            best_blocks = blocks
    return best_cut, best_blocks
```

### Reading METIS Format and Partitioning

```python
import kahip

def read_metis_graph(filename):
    """Read a METIS format graph file, return (n, xadj, adjncy, vwgt, adjcwgt)."""
    with open(filename) as f:
        lines = [l.strip() for l in f if not l.startswith('%') and l.strip()]

    header = lines[0].split()
    n, m = int(header[0]), int(header[1])
    fmt = int(header[2]) if len(header) > 2 else 0

    has_edge_weights = (fmt % 10) == 1
    has_vertex_weights = ((fmt // 10) % 10) == 1

    xadj = [0]
    adjncy = []
    adjcwgt = []
    vwgt = []

    for i, line in enumerate(lines[1:n+1]):
        tokens = list(map(int, line.split()))
        idx = 0
        if has_vertex_weights:
            vwgt.append(tokens[idx])
            idx += 1
        else:
            vwgt.append(1)

        while idx < len(tokens):
            neighbor = tokens[idx] - 1  # convert to 0-indexed
            adjncy.append(neighbor)
            idx += 1
            if has_edge_weights:
                adjcwgt.append(tokens[idx])
                idx += 1
            else:
                adjcwgt.append(1)
        xadj.append(len(adjncy))

    if not has_edge_weights:
        adjcwgt = None
    if all(w == 1 for w in vwgt):
        vwgt = None

    return n, xadj, adjncy, vwgt, adjcwgt

# Usage:
n, xadj, adjncy, vwgt, adjcwgt = read_metis_graph('examples/rgg_n_2_15_s0.graph')
edgecut, blocks = kahip.kaffpa(
    vwgt, xadj, adjcwgt, adjncy,
    nblocks=4, imbalance=0.03, suppress_output=1, seed=0, mode=kahip.ECO
)
print(f"Partitioned {n} nodes into 4 blocks, edge cut = {edgecut}")
```

### Integration with NetworkX

```python
import kahip
import networkx as nx

def partition_networkx(G, k, imbalance=0.03, mode=kahip.ECO, seed=0):
    """Partition a NetworkX graph using KaHIP."""
    nodes = list(G.nodes())
    node_to_idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)

    # Build CSR
    xadj = [0]
    adjncy = []
    adjcwgt = []
    for v in nodes:
        for u, data in G[v].items():
            adjncy.append(node_to_idx[u])
            adjcwgt.append(int(data.get('weight', 1)))
        xadj.append(len(adjncy))

    vwgt = [int(G.nodes[v].get('weight', 1)) for v in nodes]

    edgecut, blocks = kahip.kaffpa(
        vwgt, xadj, adjcwgt, adjncy, k, imbalance,
        suppress_output=1, seed=seed, mode=mode
    )

    # Map back to original node IDs
    return edgecut, {nodes[i]: blocks[i] for i in range(n)}

# Example:
G = nx.grid_2d_graph(10, 10)
edgecut, partition = partition_networkx(G, k=4)
print(f"Edge cut: {edgecut}")
print(f"Block sizes: {[sum(1 for b in partition.values() if b==i) for i in range(4)]}")
```

### Integration with SciPy Sparse Matrices

```python
import kahip
import scipy.sparse as sp
import numpy as np

def partition_scipy_sparse(A, k, imbalance=0.03, mode=kahip.ECO, seed=0):
    """Partition a graph given as a scipy sparse adjacency matrix."""
    # Convert to CSR format
    A_csr = A.tocsr()
    n = A_csr.shape[0]

    xadj   = A_csr.indptr.tolist()
    adjncy = A_csr.indices.tolist()
    adjcwgt = [int(w) for w in A_csr.data]  # edge weights from matrix values

    edgecut, blocks = kahip.kaffpa(
        None,        # uniform node weights
        xadj, adjcwgt, adjncy,
        k, imbalance, suppress_output=1, seed=seed, mode=mode
    )
    return edgecut, np.array(blocks)

# Example: Laplacian-like sparse matrix
A = sp.random(1000, 1000, density=0.01, format='csr')
A = (A + A.T)  # make symmetric (undirected)
A.setdiag(0)
A.eliminate_zeros()
# Ensure integer weights
A.data = np.ones(A.nnz, dtype=int)

edgecut, part = partition_scipy_sparse(A, k=4)
print(f"Edge cut: {edgecut}, partition shape: {part.shape}")
```

---

## Limitations

- The Python bindings expose only `kaffpa` (serial KaFFPa). ParHIP (distributed parallel) is not available from Python.
- Thread safety: do not call `kahip.kaffpa` concurrently from multiple threads.
- Large graphs: the Python interface adds overhead from list-to-array conversion. For performance-critical code with very large graphs, use the C/C++ library directly.
- The `kahip_graph` class currently does not support reading METIS format files; use the helper above.
