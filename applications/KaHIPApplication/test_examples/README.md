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

### `cube.mdpa` — 4 partitions, 3 trials per KaHIP strategy, 1991 nodal edges

| Strategy | Edge cut | Imbalance | Time (s) |
|---|---|---|---|
| fast | 423 | 1.0266 | 0.035 |
| eco | 421 | 1.0266 | 0.078 |
| strong | 439 | 1.0266 | 0.414 |
| fastsocial | 479 | 1.0266 | 0.035 |
| ecosocial | 432 | 1.0266 | 0.047 |
| strongsocial | 450 | 1.0266 | 0.297 |
| **metis** | **379** | **1.0266** | **0.018** |

Best edge cut: **METIS** (379).  Best KaHIP strategy: **eco** (421).  Fastest: **METIS** (0.018 s).

### `quads.mdpa` — 3 partitions, 3 trials per KaHIP strategy, 72 nodal edges

| Strategy | Edge cut | Imbalance | Time (s) |
|---|---|---|---|
| fast | 22 | 1.0800 | 0.002 |
| eco | 22 | 1.0800 | 0.004 |
| strong | 22 | 1.0800 | 0.033 |
| fastsocial | 22 | 1.0800 | 0.005 |
| ecosocial | 22 | 1.0800 | 0.004 |
| strongsocial | 23 | 1.0800 | 0.036 |
| **metis** | 34 | 1.0800 | 0.001 |

Best edge cut: **KaHIP fast/eco/strong/fastsocial/ecosocial** (22).  Fastest: **METIS** (0.001 s).

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
