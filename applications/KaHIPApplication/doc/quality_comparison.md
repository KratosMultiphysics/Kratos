# Partition Quality Comparison Tool

`python_scripts/compute_partition_quality.py` benchmarks all KaHIP
preconfigurations against each other and (when available) against METIS on any
Kratos `.mdpa` mesh.  It is the recommended starting point for choosing which
preconfiguration to use for a specific mesh type.

---

## Metrics

| Metric | Symbol | Optimisation direction | Description |
|--------|--------|------------------------|-------------|
| **Edge cut** | `edge_cut` | lower is better | Number of edges in the nodal-adjacency graph that cross partition boundaries.  This is the quantity KaHIP directly minimises. |
| **Load imbalance** | `imbalance` | lower is better, 1.0 = perfect | `max(block_size) / mean(block_size)`.  Values close to 1 mean that every partition owns roughly the same number of nodes. |
| **Partition time** | `time_s` | lower is better | Wall-clock seconds for `Execute()`.  Dominated by the KaHIP / METIS call; I/O is not counted. |

---

## Usage

```bash
# Default: cube.mdpa with 4 partitions, KaHIP strategies only
python applications/KaHIPApplication/python_scripts/compute_partition_quality.py

# Custom mesh and partition count
python .../compute_partition_quality.py \
    --mesh path/to/mesh.mdpa \
    --partitions 8

# Also benchmark social-graph preconfigurations
python .../compute_partition_quality.py --include-social

# More KaHIP trials per strategy (higher quality, slower)
python .../compute_partition_quality.py --num-trials 5 --partitions 8

# Adjust allowed imbalance
python .../compute_partition_quality.py --imbalance 0.05

# Custom output file
python .../compute_partition_quality.py --output results/my_mesh_quality.png
```

All options:

```
--mesh PATH       Path to the .mdpa mesh file (default: test_examples/cube.mdpa)
--partitions K    Number of partitions k (default: 4, minimum: 2)
--output FILE     Output PNG filename (default: partition_quality.png)
--imbalance F     Allowed imbalance fraction, e.g. 0.03 = 3 % (default: 0.03)
--seed INT        Random seed (default: 0)
--num-trials INT  Run KaHIP this many times per strategy, keep the best
                  result (default: 1)
--include-social  Also benchmark fastsocial / ecosocial / strongsocial
--dimension INT   Spatial dimension forwarded to METIS (2 or 3, default: 3)
```

---

## Output

The script writes a three-panel PNG figure:

```
┌─────────────────┬────────────────────────┬──────────────────┐
│   Edge cut      │   Load imbalance       │  Partition time  │
│  (bar chart)    │   (bar chart)          │  (bar chart)     │
└─────────────────┴────────────────────────┴──────────────────┘
```

It also prints a table to stdout:

```
============================================================
  KaHIP Partition Quality Comparison
============================================================
  Mesh       : .../test_examples/cube.mdpa
  Partitions : 4
  Nodal edges: 1991

Strategy        Edge cut     Imbalance     Time (s)
------------------------------------------------------
fast                 423        1.0266        0.035
eco                  421        1.0266        0.078
strong               439        1.0266        0.414
fastsocial           479        1.0266        0.035
ecosocial            432        1.0266        0.047
strongsocial         450        1.0266        0.297
metis                379        1.0266        0.018
------------------------------------------------------
  Best edge cut  : KaHIP eco   (421 cut edges)
  Best balance   : all strategies tied (imbalance = 1.0266)
  Fastest        : METIS (0.018 s)
```

> **Note:** The results above were produced on the bundled `cube.mdpa` mesh
> (413 nodes, 1191 tetrahedra, 4 partitions, 3 trials per KaHIP strategy).
> Pre-generated plots are saved alongside the mesh in `test_examples/`.

---

## KaHIP preconfigurations compared

| Strategy | Quality | Speed | Recommended when … |
|---|---|---|---|
| `fast` | ★★☆☆☆ | ★★★★★ | Time-critical, quality is secondary |
| `eco` | ★★★☆☆ | ★★★☆☆ | Default; good balance of quality and speed |
| `strong` | ★★★★★ | ★★☆☆☆ | Offline pre-processing, maximum quality |
| `fastsocial` | ★★☆☆☆ | ★★★★★ | Scale-free / power-law degree graphs |
| `ecosocial` | ★★★☆☆ | ★★★☆☆ | Social graphs, balanced mode |
| `strongsocial` | ★★★★★ | ★★☆☆☆ | Social graphs, maximum quality |
| METIS | ★★★☆☆ | ★★★★☆ | Baseline comparison |

> **Tip:** For FEM meshes, use `--include-social` and compare side-by-side —
> standard preconfigurations almost always win on structured meshes, but the
> social variants occasionally edge ahead on unstructured or adaptive meshes.

---

## METIS comparison

When `KratosMultiphysics.MetisApplication` is importable, the script
automatically adds a METIS bar to the chart.  No extra arguments are needed.
If MetisApplication is not compiled, the comparison is silently skipped.

---

## Dependencies

| Package | Required | Install |
|---|---|---|
| KratosMultiphysics + KaHIPApplication | yes | build with CMake |
| matplotlib | yes (for the plot) | `pip install matplotlib` |
| numpy | yes (for the plot) | `pip install numpy` |
| KratosMultiphysics.MetisApplication | no | build with `-DINCLUDE_METIS=ON` |

Without matplotlib/numpy the script still runs and prints the table; the PNG
is simply not written.
