---
description: "Use when building the project, running tests, or working with VS Code tasks in Kratos Multiphysics. Covers configure scripts, environment variables, and available VS Code task shortcuts."
---

# Build System and VS Code Tasks — Kratos Multiphysics

## Primary Build Entrypoints

**Never** use ad-hoc `cmake` commands as the default path. Always use the wrapper scripts:

| Platform | Template script | Personalized copy |
|----------|----------------|-------------------|
| Linux    | `scripts/standard_configure.sh` | `build/configure.sh` |
| Windows  | `scripts/standard_configure.bat` | `build/configure.bat` |

The personalized copy (`build/configure.*`) is what the VS Code `Build` task invokes.
Copy the template, set your compilers and desired applications, and use that copy locally.

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

## VS Code Tasks

Prefer existing tasks from `.vscode/tasks.json` over custom commands.

### Build

| Task | Purpose |
|------|---------|
| `Build` | Configure + compile (build type and generator selectable) |
| `MPI Build` | Build with `USE_MPI=ON` |

### Run / Test

| Task | Purpose |
|------|---------|
| `Run Tests` | Run the full Python test suite |
| `Run CurrentFile` | Run the currently open Python file with correct env |
| `Run C++ Tests` | Run all C++ GTest suites via `run_cpp_tests.py` |
| `Run C++ Test Suite` | Run a specific GTest executable (e.g. `KratosCoreTest`) |
| `Run C++ Test Suite Filtered` | Run with `--gtest_filter=*<pattern>*` |
| `Run Current Benchmark file to JSON` | Run a benchmark and write JSON output |

If a suitable task exists, use it instead of inventing custom command flows.

## Generated / Versioned Artifacts — Do Not Edit

- `build/<BuildType>/compile_commands.json` (CMake-generated)
- CMake cache files (`CMakeCache.txt`, `cmake_install.cmake`)
- Install manifests under `build/<BuildType>/`
