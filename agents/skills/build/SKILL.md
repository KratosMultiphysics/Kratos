---
name: build
description: "Use when building Kratos or running its compiled tests. The configure scripts are the required entry point for compiling this repository."
---

# Building Kratos Multiphysics

**Never** use ad-hoc `cmake` commands. Always go through the configure scripts:

| Platform | Template script | Personalized copy |
|----------|----------------|-------------------|
| Linux    | `scripts/standard_configure.sh` | `build/configure.sh` |
| Windows  | `scripts/standard_configure.bat` | `build/configure.bat` |

Copy the template once into `build/configure.*`, set compilers/applications there, then always
build through that personalized copy.

That's sufficient for almost every task. The rest of this file is only needed when a task
requires changing the compilation mode, applications built, or another configure option.

## What the Configure Scripts Do

- Set compiler paths, Python executable, and CMake generator.
- Declare the list of applications to build via `KRATOS_APPLICATIONS`.
- Invoke CMake configure + build in one call.
- Install compiled libraries and Python modules into `bin/<BuildType>/`.

## Common Build Environment Variables

| Variable | Values / Notes |
|----------|---------------|
| `KRATOS_BUILD_TYPE` | `Release`, `RelWithDebInfo`, `FullDebug`, `Custom` |
| `KRATOS_SOURCE` | Path to the repository root |
| `KRATOS_BUILD` | Build tree root (default: `<repo>/build`) |
| `KRATOS_APPLICATIONS` | Semicolon-separated list of application paths |
| `PYTHON_EXECUTABLE` | Path to the Python interpreter |
| `CMAKE_GENERATOR` | e.g. `Ninja`, `Visual Studio 17 2022` |
| `NUMBER_OF_COMPILATION_CORES` | Parallel compile jobs |
| `PYTHONPATH` | Must include `bin/<BuildType>` at runtime |
| `LD_LIBRARY_PATH` | Shared library path (Linux) — must include `bin/<BuildType>/libs` |
| `PATH` | DLL/executable path (Windows) — must include `bin/<BuildType>/libs` |

Key CMake flags passed inside the configure scripts:

| CMake option | Purpose |
|-------------|---------|
| `KRATOS_BUILD_TESTING=ON` | Compile C++ GTest binaries |
| `KRATOS_BUILD_BENCHMARK=ON` | Compile Google Benchmark binaries |
| `USE_MPI=ON` | Enable MPI-parallel builds |
| `USE_EIGEN_MKL=ON` | Link Eigen against Intel MKL |

## Generated / Versioned Artifacts — Do Not Edit

- `build/<BuildType>/compile_commands.json` (CMake-generated)
- CMake cache files (`CMakeCache.txt`, `cmake_install.cmake`)
- Install manifests under `build/<BuildType>/`
