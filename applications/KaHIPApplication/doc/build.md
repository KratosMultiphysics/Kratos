# Building KaHIP

## Requirements

### Mandatory

| Dependency | Minimum Version | Notes |
|---|---|---|
| C++ compiler | GCC 7+, Clang 6+, MSVC 2019+ | C++11 required |
| CMake | 3.10+ | |
| Make or Ninja | any | Build system |

### Optional

| Dependency | Purpose | CMake Flag |
|---|---|---|
| OpenMPI / MPICH | ParHIP distributed partitioner, kaffpaE | auto-detected |
| OpenMP | Shared-memory parallelism within serial tools | auto-detected |
| METIS + GKlib | Faster `fast_node_ordering` program | auto-detected |
| Gurobi | ILP exact solver and ILP improver | `-DUSE_ILP=ON` |
| tcmalloc | Faster heap allocation | `-DUSE_TCMALLOC=ON` |
| pybind11 | Python bindings | `-DBUILDPYTHONMODULE=ON` |
| ccache | Build caching (faster rebuilds) | auto-detected |

---

## Using the configure.sh Script

The recommended way to build is via the `build/configure.sh` script, which wraps CMake and provides sensible defaults plus dependency checking.

```bash
cd build/
./configure.sh [OPTIONS]
```

### Common Invocations

```bash
# Standard release build (auto-detects MPI, OpenMP, Metis)
./configure.sh

# Debug build, no MPI, no native CPU optimizations
./configure.sh --build-type Debug --no-mpi --no-native

# Full build: Python bindings, 64-bit edges, tcmalloc
./configure.sh --python --64bit --tcmalloc

# Portable build installed to /usr/local
./configure.sh --no-native --prefix /usr/local

# With Ninja instead of Make
./configure.sh --generator Ninja

# ILP support (requires Gurobi)
GUROBI_HOME=/opt/gurobi1100/linux64 ./configure.sh --ilp

# Deterministic ParHIP (same seed → same result)
./configure.sh --deterministic-parhip

# CI/CD: minimal, no MPI, no deploy directory
./configure.sh --no-mpi --no-native --no-deploy

# Dry run (print cmake commands without executing)
./configure.sh --python --dry-run
```

Run `./configure.sh --help` for the full option reference.

---

## Manual CMake Build

```bash
mkdir build && cd build

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

make -j$(nproc)
```

### All CMake Options

| Option | Default | Description |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release`, `Debug`, `RelWithDebInfo`, `MinSizeRel` |
| `CMAKE_INSTALL_PREFIX` | system default | Installation root |
| `NOMPI` | `OFF` | Disable all MPI-dependent targets (ParHIP, kaffpaE) |
| `PARHIP` | `ON` | Build ParHIP (requires MPI) |
| `DETERMINISTIC_PARHIP` | `OFF` | Enforce deterministic ParHIP runs |
| `BUILDPYTHONMODULE` | `OFF` | Build pybind11 Python bindings |
| `USE_ILP` | `OFF` | Build ILP improver/exact solver (requires Gurobi) |
| `USE_TCMALLOC` | `OFF` | Link against tcmalloc |
| `64BITMODE` | `OFF` | 64-bit edge indices (`kahip_idx` = `int64_t`) |
| `NONATIVEOPTIMIZATIONS` | `OFF` | Disable `-march=native` |
| `OPTIMIZED_OUTPUT` | `OFF` | Enable `KAFFPAOUTPUT` define |

### Example Configurations

```bash
# 64-bit edge support
cmake .. -DCMAKE_BUILD_TYPE=Release -D64BITMODE=ON

# No MPI (serial-only)
cmake .. -DCMAKE_BUILD_TYPE=Release -DNOMPI=ON

# With Python bindings
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILDPYTHONMODULE=ON \
    -Dpybind11_DIR=$(python3 -m pip show pybind11 | grep Location | cut -d' ' -f2)/pybind11/share/cmake/pybind11

# ILP with Gurobi
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_ILP=ON \
    -DGUROBI_HOME=/opt/gurobi1100/linux64

# Portable build (no native CPU optimizations)
cmake .. -DCMAKE_BUILD_TYPE=Release -DNONATIVEOPTIMIZATIONS=ON
```

---

## Build Targets

### Always Built

| Target | Type | Description |
|---|---|---|
| `libkaffpa` | OBJECT | Core multilevel partitioning |
| `libmapping` | OBJECT | Process mapping algorithms |
| `libspac` | OBJECT | Edge partitioning |
| `libnodeordering` | OBJECT | Node ordering / nested dissection |
| `kahip` | SHARED | Shared library interface (`libkahip.so`) |
| `kahip_static` | STATIC | Static library interface (`libkahip_static.a`) |
| `kaffpa` | Executable | Main partitioning tool |
| `global_multisection` | Executable | Hierarchical k-way partitioning |
| `evaluator` | Executable | Partition quality evaluation |
| `edge_evaluator` | Executable | Edge partition evaluation |
| `node_separator` | Executable | Vertex separator computation |
| `label_propagation` | Executable | Label propagation clustering |
| `partition_to_vertex_separator` | Executable | Convert partition → separator |
| `graphchecker` | Executable | Graph format validation |
| `edge_partitioning` | Executable | Edge partitioning (SPAC) |
| `node_ordering` | Executable | Node ordering for sparse matrices |
| `interface_test` | Executable | Library API example |

### Built with MPI (`-DNOMPI=OFF`, default)

| Target | Description |
|---|---|
| `libkaffpa_parallel` | Parallel memetic algorithm (OBJECT library) |
| `kaffpaE` | Evolutionary parallel partitioner |

### Built with MPI + ParHIP (`-DPARHIP=ON`, default when MPI present)

| Target | Description |
|---|---|
| `parhip` | Distributed parallel partitioner |
| `parhip_interface` | ParHIP shared library |
| `parhip_interface_static` | ParHIP static library |
| `toolbox` | Partition evaluation and conversion |
| `graph2binary` | Convert METIS → binary BGF format |
| `graph2binary_external` | BGF conversion (external partitioning) |
| `readbgf` | Read binary BGF format |
| `dspac` | Distributed edge partitioning |

### Built conditionally

| Target | Condition | Description |
|---|---|---|
| `fast_node_ordering` | METIS found | Fast ordering using reductions + METIS |
| `ilp_improve` | `-DUSE_ILP=ON` | ILP-based partition improvement |
| `ilp_exact` | `-DUSE_ILP=ON` | ILP exact solver |
| `kahip_python_binding` | `-DBUILDPYTHONMODULE=ON` | Python module |

---

## Compiler Flags Applied

The build system automatically applies the following compiler flags when supported:

| Flag | Purpose |
|---|---|
| `-O3` (Release) | Full optimization |
| `-march=native` | CPU-specific optimization (disable with `NONATIVEOPTIMIZATIONS`) |
| `-funroll-loops` | Loop unrolling |
| `-fno-stack-limit` | Remove stack size limit |
| `-fpermissive` | Relaxed type checking |
| `-Wall` | All warnings |
| `-Wno-unused-result` | Suppress unused-result warnings |
| `-Wno-sign-compare` | Suppress sign-comparison warnings |

---

## Installation

After building, install to the system (or a custom prefix):

```bash
# Install to default prefix (/usr/local on Linux/macOS)
cmake --install build/

# Install to custom prefix
cmake --install build/ --prefix /opt/kahip

# Installed layout:
# /opt/kahip/
# ├── bin/          executables (kaffpa, parhip, …)
# ├── lib/          libkahip.so, libkahip_static.a, libparhip.so, …
# ├── include/      kaHIP_interface.h, parhip_interface.h
# └── lib/pkgconfig/ kahip.pc, parhip.pc
```

### Using pkg-config after installation

```bash
# Compile flags
pkg-config --cflags kahip

# Link flags
pkg-config --libs kahip

# In your build
g++ myapp.cpp $(pkg-config --cflags --libs kahip) -o myapp
```

---

## Platform Notes

### Linux

Standard build. OpenMPI or MPICH must be installed separately.

```bash
# Ubuntu/Debian
sudo apt-get install cmake libopenmpi-dev openmpi-bin

# Fedora/RHEL
sudo dnf install cmake openmpi-devel
```

### macOS

```bash
# Homebrew (recommended — installs pre-built binaries)
brew install KaHIP/kahip/kahip

# Or build from source
brew install cmake open-mpi
cd build && ./configure.sh
```

### Windows

KaHIP builds on Windows with MSVC 2019+. MPI support requires Microsoft MPI (MS-MPI).

```powershell
# Using cmake directly
cmake .. -DCMAKE_BUILD_TYPE=Release -DNOMPI=ON
cmake --build . --config Release
```

---

## Troubleshooting

### MPI not found

```
CMake Error: Could not find MPI.
```

Set `MPI_HOME` or use `--no-mpi`:

```bash
MPI_HOME=/usr/lib/x86_64-linux-gnu/openmpi cmake ..
# or
cmake .. -DNOMPI=ON
```

### pybind11 not found for Python bindings

```bash
pip install pybind11
export pybind11_DIR=$(python3 -m pip show pybind11 | grep Location | cut -d' ' -f2)/pybind11/share/cmake/pybind11
cmake .. -DBUILDPYTHONMODULE=ON
```

### Gurobi not found for ILP

```bash
export GUROBI_HOME=/opt/gurobi1100/linux64
cmake .. -DUSE_ILP=ON -DGUROBI_HOME=$GUROBI_HOME
```

### "undefined reference to tcmalloc"

tcmalloc is not installed. Either install `google-perftools` or remove `-DUSE_TCMALLOC=ON`.

### Stack overflow at runtime

KaHIP uses deep recursion. If you get stack overflows on large graphs, the `-fno-stack-limit` flag should help. On Linux you can also increase the stack size:

```bash
ulimit -s unlimited
```
