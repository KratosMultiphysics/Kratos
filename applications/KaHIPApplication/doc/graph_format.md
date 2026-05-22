# Graph File Formats

KaHIP uses two graph file formats: the **METIS text format** (default, human-readable) and the **binary BGF format** (fast parallel I/O for ParHIP).

---

## METIS Text Format

File extension: `.graph`

### Structure

```
<n> <m> [<fmt>] [<ncon>]
<v1_adj1> [<v1_adj2> ...]
<v2_adj1> [<v2_adj2> ...]
...
```

Where:
- `n` = number of vertices
- `m` = number of edges (undirected, counted once)
- `fmt` = format specifier (optional, see below)
- `ncon` = number of vertex weights per node (optional, default 1 if vertex weights present)
- Each subsequent line describes one vertex's adjacency list (1-indexed, space-separated)

Lines starting with `%` are comments.

### Format Specifier

The `fmt` field is a 3-digit number `xyz` where:
- `z` = 0: no edge weights; 1: edge weights present
- `y` = 0: no vertex weights; 1: vertex weights present
- `x` = 0: no vertex sizes; 1: vertex sizes present

| `fmt` value | Meaning |
|---|---|
| `0` or absent | Unweighted graph |
| `1` | Edge-weighted graph |
| `10` | Vertex-weighted graph |
| `11` | Both vertex and edge weights |
| `100` | Vertex sizes only |

### Unweighted Graph Example

5-node graph with edges: 0-1, 0-4, 1-2, 1-4, 2-3, 3-4

```
5 6
2 5        ← node 1 (1-indexed) connects to nodes 2 and 5
1 3 5      ← node 2 connects to nodes 1, 3, 5
2 4        ← node 3 connects to nodes 2 and 4
3 5        ← node 4 connects to nodes 3 and 5
1 2 4      ← node 5 connects to nodes 1, 2, 4
```

Note: METIS format uses **1-indexed** node IDs. KaHIP internally converts to 0-indexed.

### Edge-Weighted Graph Example (`fmt = 1`)

```
5 6 1
2 1 5 1          ← node 1: neighbor 2 (weight 1), neighbor 5 (weight 1)
1 1 3 2 5 1      ← node 2: neighbor 1 (w=1), neighbor 3 (w=2), neighbor 5 (w=1)
2 2 4 3          ← node 3: neighbor 2 (w=2), neighbor 4 (w=3)
3 3 5 1          ← node 4: neighbor 3 (w=3), neighbor 5 (w=1)
1 1 2 1 4 1      ← node 5: neighbor 1 (w=1), neighbor 2 (w=1), neighbor 4 (w=1)
```

Each neighbor entry is `<node_id> <edge_weight>`.

### Vertex-Weighted Graph Example (`fmt = 10`)

```
5 6 10
3 2 5              ← node 1 has weight 3, connects to 2 and 5
1 1 3 5            ← node 2 has weight 1, connects to 1, 3, 5
2 2 4              ← node 3 has weight 2, connects to 2 and 4
1 3 5              ← node 4 has weight 1, connects to 3 and 5
2 1 2 4            ← node 5 has weight 2, connects to 1, 2, 4
```

Format per line: `<vertex_weight> <neighbor1> <neighbor2> ...`

### Both Weighted (`fmt = 11`)

```
5 6 11
3 2 1 5 1          ← node 1: vwgt=3, connects to 2 (ew=1), 5 (ew=1)
1 1 1 3 2 5 1      ← node 2: vwgt=1, neighbor 1 (ew=1), 3 (ew=2), 5 (ew=1)
2 2 2 4 3          ← node 3: vwgt=2, neighbor 2 (ew=2), 4 (ew=3)
1 3 3 5 1          ← node 4: vwgt=1, neighbor 3 (ew=3), 5 (ew=1)
2 1 1 2 1 4 1      ← node 5: vwgt=2, neighbor 1 (ew=1), 2 (ew=1), 4 (ew=1)
```

### Fast I/O (`--mmap_io`)

For large graphs, the `kaffpa` executable supports `--mmap_io` which uses memory-mapped file access instead of `fscanf`. This can be an order of magnitude faster.

---

## Binary BGF Format

File extension: `.bgf` (Binary Graph Format)

BGF is optimized for fast parallel reading in ParHIP. It stores the graph in a compact binary layout with metadata for efficient distributed loading.

### Converting to BGF

```bash
# Convert METIS text format to binary BGF
./graph2binary input.graph output.bgf

# Convert with an existing partition
./graph2binary_external input.graph partition.txt output.bgf
```

### When to Use BGF

- Always use BGF when running ParHIP on large graphs (millions of nodes)
- BGF enables truly parallel I/O: each rank reads only its portion
- Text format requires sequential parsing, which bottlenecks on many ranks

---

## Partition File Format

Partition output files (e.g., written with `--output_filename`) use a simple text format:

```
<block_of_node_0>
<block_of_node_1>
...
<block_of_node_n-1>
```

One integer per line, 0-indexed block IDs. Example partition of 5 nodes into 2 blocks:

```
0
0
1
1
0
```

---

## Graph Validation

Use `graphchecker` to validate a graph file:

```bash
./graphchecker graph.graph
```

This checks:
- Correct header format
- Consistent edge count
- All edges appear in both directions (undirected requirement)
- Edge weights (if present) are consistent
- Self-loops and multi-edges (these are generally not supported)

---

## CSR Format (Library Interface)

When using the C/C++ API directly, the graph is passed as CSR arrays (see [api_serial.md](api_serial.md)):

```
xadj[0..n]        row offsets: node i's neighbors are adjncy[xadj[i]..xadj[i+1]-1]
adjncy[0..2m-1]   neighbor node IDs (0-indexed)
vwgt[0..n-1]      node weights (NULL = uniform 1)
adjcwgt[0..2m-1]  edge weights (NULL = uniform 1)
```

**Conversion from METIS format to CSR**:

```cpp
// METIS is 1-indexed; CSR is 0-indexed
// METIS line i: "w_i adj1 adj2 ..." → xadj[i-1]..xadj[i], adjncy entries
```

### CSR Example (same 5-node graph, 0-indexed)

```
n = 5, m = 6 (12 directed edge entries)

xadj   = [0, 2, 5, 7, 9, 12]
adjncy = [1, 4,      ← node 0 connects to 1, 4
          0, 2, 4,   ← node 1 connects to 0, 2, 4
          1, 3,      ← node 2 connects to 1, 3
          2, 4,      ← node 3 connects to 2, 4
          0, 1, 3]   ← node 4 connects to 0, 1, 3
```

---

## Example Graphs

The repository includes several example graphs in `examples/`:

| File | Nodes | Edges | Description |
|---|---|---|---|
| `rgg_n_2_15_s0.graph` | 32,768 | ~160K | Random geometric graph (mesh-like) |
| `rgg_n_2_15_s0.bgf` | 32,768 | ~160K | Same graph in binary BGF format |
| `delaunay_n15.graph` | 32,768 | ~98K | Delaunay triangulation graph |
| `example_weighted.graph` | small | small | Small graph with both vertex and edge weights |

These are suitable for quick testing:
```bash
./kaffpa examples/rgg_n_2_15_s0.graph --k 4 --preconfiguration=strong
```

---

## Common Graph Sources

| Source | URL / Tool |
|---|---|
| SuiteSparse Matrix Collection | https://sparse.tamu.edu/ (download in METIS format) |
| SNAP (Stanford) | http://snap.stanford.edu/data/ (need conversion) |
| 10th DIMACS Challenge | http://www.cc.gatech.edu/dimacs10/ |
| KaGen (graph generator) | https://github.com/KaHIP/KaGen |

### Converting SNAP Format to METIS

SNAP uses a simple edge list format. Convert with:
```bash
# One-liner: convert edge list to METIS format
# (assuming 0-indexed nodes in edge list)
awk 'NF==2 {print $1+1, $2+1}' snap_graph.txt | \
    awk 'BEGIN{n=0;m=0}{if($1>n)n=$1;if($2>n)n=$2;edges[m++]=$0}
         END{print n,m/2; for(i=0;i<m;i++) print edges[i]}' > metis_graph.txt
```

For reliable conversion, use KaGen or a dedicated converter.
