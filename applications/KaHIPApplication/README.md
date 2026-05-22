# KaHIPApplication

KaHIP-based graph-partitioning application for Kratos Multiphysics.

## What it does

`KaHIPApplication` partitions a Kratos `ModelPart` (mesh) for parallel computation,
serving as a drop-in alternative to `MetisApplication`. It:

1. Reads the nodal connectivity graph from the mesh via the `IO` interface.
2. Converts the graph to CSR format.
3. Calls [KaHIP](https://github.com/KaHIP/KaHIP) (`kaffpa()`) to partition the nodal
   graph into *k* blocks.
4. Assigns elements, conditions, geometries, and master-slave constraints to partitions
   based on the majority-partition of their nodes (with optional boundary synchronisation).
5. Writes per-rank `.mdpa` files via `ModelPartIO::DivideInputToPartitions()`.

## Why KaHIP instead of METIS?

- KaHIP typically produces lower edge cuts (fewer inter-partition communication
  edges), especially with the `eco` and `strong` preconfigurations.
- `eco` is comparable in speed to METIS's k-way partitioner while producing
  5–20 % fewer edge cuts on FEM meshes.
- `strong` can reduce edge cuts further at the cost of longer runtime.
- KaHIP is MIT licensed.

## Prerequisites

- Kratos Multiphysics (built with this application enabled)
- KaHIP ≥ v3.25 (the in-tree build under `external_libraries/KaHIP` is used by default)
- OpenMP (required by KaHIP)

## CMake Configuration

```bash
# Enable the application (add to configure.sh or cmake command line)
add_app ${KRATOS_APP_DIR}/KaHIPApplication

# Optional: point to a system or custom KaHIP installation
-DKAHIP_ROOT=/path/to/kahip/build

# Optional: enable 64-bit edge indices (for meshes with > 2^31 edge entries)
-DKAHIP_64BIT=ON

# Optional: disable ParHIP (MPI-parallel partitioning)
-DKAHIP_USE_PARHIP=OFF
```

## Python API

### Drop-in replacement for MetisApplication

```python
import KratosMultiphysics as KM
import KratosMultiphysics.KaHIPApplication as KaHIP

model = KM.Model()
io_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER
io = KM.ModelPartIO("my_mesh", io_flags)

# Same positional arguments as MetisDivideHeterogeneousInputProcess
partitioner = KaHIP.KaHIPDivideHeterogeneousInputProcess(
    io,
    number_of_partitions=4,
    dimension=3,
    verbosity=0,
    synchronize_conditions=True)

partitioner.Execute()
# Writes my_mesh_partitioned/my_mesh_0.mdpa … my_mesh_3.mdpa
```

### Using Parameters for fine-grained control

```python
import KratosMultiphysics as KM
import KratosMultiphysics.KaHIPApplication as KaHIP

io = KM.ModelPartIO("my_mesh", KM.ModelPartIO.READ)

settings = KM.Parameters("""{
    "preconfiguration":     "eco",
    "imbalance":             0.03,
    "seed":                  0,
    "suppress_output":       true,
    "num_trials":            3,
    "verbosity":             1,
    "synchronize_conditions": true
}""")

partitioner = KaHIP.KaHIPDivideHeterogeneousInputProcess(io, 4, settings)
partitioner.Execute()
```

## Partitioner Settings

| Key                      | Type    | Default | Description                                                  |
|--------------------------|---------|---------|--------------------------------------------------------------|
| `preconfiguration`       | string  | `"eco"` | KaHIP mode: `fast`, `eco`, `strong`, `fastsocial`, `ecosocial`, `strongsocial` |
| `imbalance`              | double  | `0.03`  | Allowed weight imbalance fraction (0.03 = 3 %)             |
| `seed`                   | int     | `0`     | Base random seed; trial i uses seed+i                       |
| `suppress_output`        | bool    | `true`  | Suppress KaHIP's stdout                                      |
| `num_trials`             | int     | `1`     | Run KaHIP this many times; keep the best (lowest edge cut)  |
| `verbosity`              | int     | `0`     | 0=silent, 1=info, 2=debug, 3=full debug                     |
| `synchronize_conditions` | bool    | `false` | Co-locate boundary conditions with their parent element      |

## Running the tests

```bash
# Python tests
cd path/to/kratos/build
python -m pytest applications/KaHIPApplication/tests/

# C++ unit tests (requires KRATOS_BUILD_TESTING=ON)
ctest -R KaHIPApplication --verbose
```
