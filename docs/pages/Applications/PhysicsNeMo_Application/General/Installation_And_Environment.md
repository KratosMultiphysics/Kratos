---
title: Installation and environment
keywords: installation build configure pip torch physicsnemo cuda onnxruntime cupy mkl pep 668 abi stale
tags: [Installation_And_Environment.md]
sidebar: physicsnemo_application
summary: Building the application, the runtime environment, every optional dependency with the version verified here, the CUDA-specific traps, and the one Kratos build trap that silently breaks unrelated applications.
---

# Installation and environment

## Building the application

The application compiles like any other Kratos application - its C++ side registers nothing and links only against `KratosCore`:

```bash
# in your copy of scripts/standard_configure.sh
add_app ${KRATOS_APP_DIR}/PhysicsNeMoApplication
```

There is deliberately **no CMake option** for torch or physicsnemo, no `find_package(Torch)`, no CUDA detection. The ML stack is a pure Python runtime dependency resolved lazily by every module; adding such a gate is the one contribution the [Overview](Overview.html) asks you not to make.

Several features run real Kratos solves and therefore need other applications compiled. What each needs:

| Feature | Kratos applications |
|---|---|
| the thermal cases (most tests, notebooks 07 and later, the Examples) | `ConvectionDiffusionApplication`, `LinearSolversApplication` |
| structural and transient structural cases, the adjoint cross-validation | `StructuralMechanicsApplication` (+ `ConstitutiveLawsApplication` for the thermo-mechanical case) |
| the lid-driven cavity | `FluidDynamicsApplication` |
| ROM surrogates | `RomApplication` (to *produce* a basis; consuming one needs nothing compiled) |
| adaptive remeshing | `MeshingApplication` with MMG |
| the co-simulation surrogate | `CoSimulationApplication`, `MappingApplication` |
| non-matching field transfer | `MappingApplication` |
| mesh smoothing in shape optimization | `MeshMovingApplication` |
| NURBS modeler workflows | `IgaApplication` (the sampling itself needs only the core) |
| MPI tests | `MetisApplication`, `TrilinosApplication`, an MPI build |

Every such test is conditioned on `CheckIfApplicationsAvailable` and skips with a message naming what is missing.

## Runtime environment

```bash
export PYTHONPATH=/path/to/Kratos/bin/Release
export LD_LIBRARY_PATH=/path/to/Kratos/bin/Release/libs:$LD_LIBRARY_PATH
# if StructuralMechanics' adjoint stack (MKL-linked Eigen solvers) is used:
export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib:/opt/intel/oneapi/compiler/latest/lib:$LD_LIBRARY_PATH
```

A missing MKL path shows up as `libmkl_rt.so.2: cannot open shared object file` deep inside a solve; it looks like a build failure and is not one.

## The Python dependencies

None is required for `import KratosMultiphysics.PhysicsNeMoApplication` to succeed. Each is checked lazily at the entry point that needs it, and the error names the package to install.

```bash
pip install torch                 # CPU or CUDA build, your choice; match your driver
pip install nvidia-physicsnemo    # 2.2.0 is what this application tracks
```

| Package | Verified here | Unlocks | Notes |
|---|---|---|---|
| `torch` | 2.13.0+cu130 | everything ML | |
| `nvidia-physicsnemo` | 2.2.0 | everything ML | see [Versions and compatibility](../PhysicsNeMo_Basics/Versions_And_Compatibility.html) for the extras |
| `torch_geometric`, `torch_scatter` | 2.8.0, 2.1.2 | the graph processes | `torch_scatter` is a slow source build; upstream checks both |
| `torch_sparse` or `dgl` | absent | GraphCast | no torch-2.13 wheel of `torch_sparse` was available; `dgl` predates torch 2.13 |
| `onnxscript`, `onnxruntime` | 0.7.1, 1.29.0 | ONNX export and CPU inference | torch's exporter needs `onnxscript` |
| `onnxruntime-gpu` | 1.29.0 (in a throwaway venv) | ONNX on the GPU | **replaces** `onnxruntime` - both install the same package directory and overwrite each other file for file; keep one per environment |
| `tritonclient[http]` or `[grpc]` | 2.71.0 | Triton serving | without a protocol extra it raises `RuntimeError`, not `ImportError` |
| `gpytorch` | 1.15.2 | the GP uncertainty head | or `nvidia-physicsnemo[uq-extras]` |
| `pyvista` | 0.48.4 | mesh rendering, `cfd_bridge`, the notebook figures | |
| `nvidia-physicsnemo-cfd` | 0.0.3a0 | CFD metrics, hybrid initialization | **not on PyPI**: `pip install git+https://github.com/NVIDIA/physicsnemo-cfd` |
| `cupy-cuda13x` (or `12x`) | 14.1.1 | the opt-in array backend | must match the CUDA toolkit torch was built for |
| `usd-core` | 26.8 | OpenUSD export | no Omniverse install needed to *write* a stage |
| `tetgen` | 0.8.4 | exact boundary recovery in the tetrahedral fill | AGPL; only used when `"method": "tetgen"` is asked for |
| `scipy` | | periodic particle graphs, the bistride hierarchy, the adjoint solves | pulled in by most of the above |
| `matplotlib`, `jupyterlab` | | the notebooks | |
| `physicsnemo-curator` | declined | the curator sinks | git-only, and its build downloads a Rust toolchain; the roadmap's curator-free Zarr export is the way around it |

On a PEP 668 "externally managed" interpreter (Debian and Ubuntu system Python), `pip install` needs `--break-system-packages` or a virtual environment. If you use a venv, add a `.pth` file pointing at the site-packages that hold torch and physicsnemo, or install them into the venv.

## CUDA notes

- `"device": "auto"` picks CUDA when torch sees a GPU. Every process accepts `"cpu"`, `"cuda"` or `"cuda:N"`.
- The Warp backend of several upstream functionals (SDF, deformers) is auto-selected on a CUDA machine and computes in **float32**; the bridges pin the torch backend for float64 Kratos fields where it matters.
- ONNX Runtime reports two failures as success - a silent fallback to CPU and a dropped device index. `"require_device": true` turns the first into an error; the second is fixed in the bridge.
- Two ranks cannot share one GPU under NCCL; the MPI tests run over gloo on CPU and are the reason the multi-GPU transport is a roadmap item.
- `DistributedManager.initialize_mesh` and `fully_shard()` build a **CUDA** device mesh whenever a GPU is visible, even for a gloo run; the distributed utilities here build the CPU mesh explicitly.

## Verifying an installation

```bash
cd applications/PhysicsNeMoApplication/tests
python3 test_import_contract.py          # the application imports with none of the ML packages
python3 test_suite_registration.py       # documentation and registration guards, torch-free
python3 test_PhysicsNeMoApplication.py -l small   # the serial suite; read the skip reasons
OMP_NUM_THREADS=1 mpiexec -np 2 python3 test_PhysicsNeMoApplication_mpi.py --using-mpi
```

Where this was written, the serial suite is 882 tests with 22 skips (all of them honest optional-dependency gates) and the MPI suite 37 tests at two and three ranks. If your skip count is higher, the skip messages say which package or application is missing.

## A Kratos build trap worth knowing before you add this application

Kratos reads `KRATOS_APPLICATIONS` from the *environment* and never writes it to the CMake cache, and `make install` never prunes the install tree. An install directory therefore accumulates every application ever built, while the build graph covers only the applications in the *current* configure script - and the two drift apart silently. Applications outside the graph are not rebuilt, so the next relink of `libKratosCore.so` leaves them ABI-stale: they still exist, `CheckIfApplicationsAvailable` still reports them present, and importing them fails with an undefined symbol or a segfault. Because several Kratos applications import their optional companions on a *presence* check rather than a working one, one stale library can take down an application that was never touched.

So: reconfiguring to add an application is not a local operation. Restore the **full** historical application list in the same run, and treat the install tree as unreproducible from the cache alone - nothing records what produced it. A clean rebuild takes about four hours on a laptop; the three symptoms above are cheaper to recognise than to debug.

Next: [Testing and contributing](Testing_And_Contributing.html), or [From scratch](From_Scratch.html).
