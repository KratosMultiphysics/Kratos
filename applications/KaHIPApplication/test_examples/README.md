# KaHIPApplication — Test Examples

Example meshes for testing and benchmarking the KaHIPApplication.

## Files

| File | Nodes | Elements | Conditions | Description |
|---|---|---|---|---|
| `cube.mdpa` | 413 | 1191 | 780 | 3-D tetrahedral cube mesh |
| `quads.mdpa` | 25 | 16 | — | 2-D quadrilateral mesh (quick tests) |
| `cube_partition_quality.png` | — | — | — | Pre-generated quality comparison plot (4 partitions) |
| `quads_partition_quality.png` | — | — | — | Pre-generated quality comparison plot (3 partitions) |

## Pre-generated benchmark results

### `cube.mdpa` — 8 partitions, 3 trials per KaHIP strategy, 1991 nodal edges

| Strategy | Edge cut | Imbalance | Time (s) |
|---|---|---|---|
| fast | 646 | 1.0266 | 0.057 |
| eco | 642 | 1.0266 | 0.148 |
| **strong** | **627** | **1.0266** | 0.808 |
| fastsocial | 714 | 1.0266 | 0.057 |
| ecosocial | 718 | 1.0073 | 0.092 |
| strongsocial | 664 | 1.0266 | 0.368 |
| metis | 653 | 1.0460 | 0.024 |

Best edge cut: **KaHIP Strong** (627).  Best balance: **KaHIP EcoSocial** (1.0073).  Fastest: **METIS** (0.024 s).

### `quads.mdpa` — 8 partitions, 3 trials per KaHIP strategy, 72 nodal edges

| Strategy | Edge cut | Imbalance | Time (s) |
|---|---|---|---|
| fast | 44 | 1.2800 | 0.003 |
| **eco** | **43** | **1.2800** | 0.008 |
| **strong** | **43** | **1.2800** | 0.032 |
| fastsocial | 44 | 1.2800 | 0.008 |
| **ecosocial** | **43** | **1.2800** | 0.009 |
| **strongsocial** | **43** | **1.2800** | 0.032 |
| metis | 44 | 1.9200 | 0.001 |

Best edge cut: **KaHIP eco/strong/ecosocial/strongsocial** (43).  Best balance: **KaHIP** strategies (1.2800 vs METIS 1.9200).  Fastest: **METIS** (0.001 s).

> Results produced on the bundled meshes with all 6 KaHIP preconfigurations
> and MetisApplication as a baseline.  See `cube_partition_quality.png` and
> `quads_partition_quality.png` for the three-panel bar charts.

## Running the benchmark yourself

```bash
# Default: cube.mdpa, 4 partitions, KaHIP strategies only
python applications/KaHIPApplication/python_scripts/compute_partition_quality.py

# Custom mesh and partition count
python .../compute_partition_quality.py --mesh path/to/my_mesh.mdpa --partitions 8

# Include social-graph preconfigurations and run 3 trials
python .../compute_partition_quality.py --include-social --num-trials 3

# Save figure to a custom path
python .../compute_partition_quality.py --output results/quality.png
```

See [`doc/quality_comparison.md`](../doc/quality_comparison.md) for the full
option reference.

## Notes

- These meshes are identical to the ones in `MetisApplication/tests/`.
- For production benchmarks, use larger meshes (tens of thousands of elements).
