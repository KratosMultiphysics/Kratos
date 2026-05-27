# KaHIPApplication — Test Examples

Example meshes for testing and benchmarking the KaHIPApplication.

## Files

| File | Description |
|---|---|
| `cube.mdpa` | 3-D tetrahedral cube mesh (413 nodes, 1191 elements, 780 conditions) |
| `quads.mdpa` | 2-D quadrilateral mesh (small, useful for quick tests) |

## Partition quality comparison

The `python_scripts/compute_partition_quality.py` utility benchmarks all KaHIP
preconfigurations and (optionally) METIS on any `.mdpa` mesh.  It uses
`cube.mdpa` by default:

```bash
# Serial run — KaHIP strategies only (default 4 partitions)
python applications/KaHIPApplication/python_scripts/compute_partition_quality.py

# Custom mesh and partition count
python .../compute_partition_quality.py --mesh path/to/my_mesh.mdpa --partitions 8

# Include social-graph preconfigurations
python .../compute_partition_quality.py --include-social

# More trials for higher quality (slower)
python .../compute_partition_quality.py --num-trials 5 --partitions 8

# METIS comparison (requires MetisApplication to be compiled)
python .../compute_partition_quality.py
```

The script writes a `partition_quality.png` figure to the current working
directory comparing edge cut, load imbalance, and partition time across
strategies.

## Notes

- These meshes are identical to the ones in `MetisApplication/tests/`.
- The `cube.mdpa` mesh is a good representative 3-D FEM mesh for comparing
  partitioning quality between strategies.
- For production benchmarks, use larger meshes (tens of thousands of elements).
