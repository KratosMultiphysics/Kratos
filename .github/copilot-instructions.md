# GitHub Copilot Instructions — Kratos Multiphysics

## 📌 Project Overview

**Kratos Multiphysics** is an open-source framework for building parallel,
multi-disciplinary simulation software. It provides a core finite-element
engine (`kratos/`) plus a library of domain-specific **applications**
(structural mechanics, fluid dynamics, contact, etc.) and exposes its
full API to Python via **pybind11**.

---

### 🗂️ Repository Structure

```
.
├── .github/
│   ├── copilot-instructions.md        # This file
│   └── workflows/                     # GitHub Actions CI definitions
│       ├── ci.yml
│       └── nightly_build.yml
├── kratos/                            # Core framework (C++ & Python)
│   ├── includes/                      # Public headers (checks.h, expect.h, define.h, …)
│   ├── sources/                       # Core C++ implementation
│   ├── containers/
│   ├── elements/
│   ├── conditions/
│   ├── processes/
│   ├── utilities/
│   ├── linear_solvers/
│   ├── solving_strategies/
│   ├── geometries/
│   ├── integration/
│   ├── spatial_containers/
│   ├── input_output/
│   ├── python/                        # pybind11 bindings for core
│   ├── python_scripts/                # Core Python modules
│   ├── testing/                       # Testing framework (testing.h, …)
│   ├── tests/                         # Core C++ and Python tests
│   ├── mpi/                           # MPI-aware core components
│   ├── benchmarks/                    # Core Google Benchmark files
│   └── CMakeLists.txt
├── applications/                      # Domain-specific Kratos applications
│   └── <ApplicationName>/
│       ├── <application_name>_application.h
│       ├── <application_name>_application.cpp
│       ├── <application_name>_application_variables.h/.cpp
│       ├── <ApplicationName>.py       # Python __init__
│       ├── CMakeLists.txt
│       ├── custom_elements/
│       ├── custom_conditions/
│       ├── custom_constitutive/
│       ├── custom_processes/
│       ├── custom_utilities/
│       ├── custom_strategies/
│       ├── custom_python/
│       │   ├── <app>_python_application.cpp  # PYBIND11_MODULE
│       │   └── add_custom_*_to_python.cpp    # Per-category bindings
│       ├── python_scripts/
│       └── tests/
│           ├── test_<ApplicationName>.py     # Suite entry point
│           └── cpp_tests/
│               ├── <app>_fast_suite.h/.cpp   # GTest fixture
│               └── test_*.cpp
├── external_libraries/                # Third-party dependencies (amgcl, …)
├── scripts/                           # Configure scripts for various platforms
│   ├── standard_configure.sh          # Linux configure + build
│   ├── standard_configure.bat         # Windows configure + build
│   └── ...
├── cmake_modules/                     # Custom CMake Find/utility modules
├── CMakeLists.txt                     # Root CMake entry point
├── INSTALL.md                         # Build instructions
├── CONTRIBUTING.md                    # Contribution guidelines
└── ...
```

---

### 🛠️ Tech Stack

| Layer            | Technology                                                    |
|------------------|---------------------------------------------------------------|
| Core language    | C++20                                                         |
| Scripting        | Python 3.x                                                    |
| Python bindings  | pybind11                                                      |
| Build system     | CMake (wrapped by `scripts/standard_configure.*`)             |
| Testing (C++)    | Google Test (GTest) + Kratos macros (`KRATOS_EXPECT_*`)       |
| Testing (Python) | `KratosUnittest` (wrapper over `unittest`)                    |
| Benchmarking     | Google Benchmark                                              |
| CI/CD            | GitHub Actions (`.github/workflows/`)                         |
| Dependencies     | Bundled in `external_libraries/`                              |

---

### ⚙️ Build System

- **Linux**: Copy and customize `scripts/standard_configure.sh`.
- **Windows**: Copy and customize `scripts/standard_configure.bat`.
- Both scripts wrap **CMake** — they set compilers, select applications,
  configure, and build in one invocation.
- The local workspace uses **`build/configure.bat`** (Windows) or
  **`build/configure.sh`** (Linux) as a personalized configure script
  invoked by the VS Code `Build` task.
- Key CMake options:
  - `KRATOS_BUILD_TYPE` — `Release` | `RelWithDebInfo` | `FullDebug` | `Custom`
  - `KRATOS_BUILD_TESTING=ON` — compile C++ GTest binaries
  - `KRATOS_BUILD_BENCHMARK=ON` — compile Google Benchmark binaries
  - `USE_MPI=ON` — enable MPI-parallel builds
  - `USE_EIGEN_MKL=ON` — link Eigen against Intel MKL
  - `KRATOS_APPLICATIONS` — semicolon-separated list of application paths

---

### 🧱 Architecture & Design Principles

#### Core vs Applications
- **`kratos/`** contains the framework core: containers, solvers, geometries,
  I/O, utilities, and the Python binding layer.
- **`applications/`** extend the core with domain-specific Elements, Conditions,
  ConstitutiveLaws, Processes, and strategies.
- Each application registers its components in `Register()` /
  `RegisterKratosApplication()` methods.

#### C++ Conventions
- Follow **Kratos coding style** for all C++ code:
  - Use `PascalCase` for class names (e.g., `TotalLagrangian`).
  - Use `PascalCase` for methods (e.g., `CalculateLocalSystem`, `GetDofList`).
  - Use `snake_case` for file names (e.g., `total_lagrangian.h`).
  - Prefix member variables with `m` (e.g., `mThickness`, `mConstitutiveLaw`).
  - Prefix reference arguments with `r` (e.g., `rModelPart`, `rCurrentProcessInfo`).
  - Prefix pointer arguments with `p` (e.g., `pElement`, `pNode`).
  - Use `const` wherever possible.
  - Prefer smart pointers (`Kratos::shared_ptr`, `Kratos::unique_ptr`) over raw pointers.
  - Use Kratos type aliases (`IndexType`, `SizeType`, `MatrixType`, `VectorType`, etc.)
    instead of raw STL or primitive types where applicable.

#### Python Conventions
- Python scripts follow **PEP 8**.
- Use `snake_case` for all Python identifiers.
- Import from compiled modules via:
  `import KratosMultiphysics` and
  `import KratosMultiphysics.<ApplicationName>Application`.
- Analysis stages and processes should inherit from the appropriate
  Kratos base classes.

#### pybind11 Bindings
- All C++ ↔ Python bindings are defined in `custom_python/` directories within
  each application (or `kratos/python/` for core).
- Binding files are named following the pattern `add_custom_<type>_to_python.cpp`
  (e.g., `add_custom_processes_to_python.cpp`).
- Use `py::class_<Derived, Derived::Pointer, Base>` with appropriate base classes.
- Always provide docstrings for exposed classes and methods.
- Use `py::return_value_policy` explicitly when returning references or pointers.

---

### 🧪 Testing

#### Python Tests (KratosUnittest)
- All Python tests use **`KratosMultiphysics.KratosUnittest`**, a wrapper
  around `unittest` with Kratos-specific helpers.
- Test files are located in `tests/` directories within each application
  (core tests in `kratos/tests/`).
- Test files are named `test_<feature_name>.py`.
- Test classes inherit from `KratosUnittest.TestCase`, which provides:
  - `skipTestIfApplicationsNotAvailable(*apps)`
  - `assertVectorAlmostEqual(v1, v2, places=7)`
  - `assertMatrixAlmostEqual(m1, m2, places=7)`
- Each application has a suite entry point
  `test_<ApplicationName>.py` that assembles suites:
  `small`, `nightly`, `validation`, `all` (plus `mpi_*` variants).
- Test factory patterns are common — a base class runs an analysis
  from a JSON `Parameters` file; concrete tests only set `file_name`.

  ```python
  import KratosMultiphysics
  import KratosMultiphysics.KratosUnittest as KratosUnittest

  class TestMyFeature(KratosUnittest.TestCase):

      def setUp(self):
          self.model = KratosMultiphysics.Model()
          self.model_part = self.model.CreateModelPart("TestPart")

      def test_something(self):
          # Arrange / Act / Assert
          self.assertAlmostEqual(expected, actual, places=6)

  if __name__ == "__main__":
      KratosUnittest.main()
  ```

#### C++ Tests (GTest + Kratos Wrappers)
- All C++ tests use **Google Test** via Kratos wrapper macros from
  `testing/testing.h`.
- The primary test macro is **`KRATOS_TEST_CASE_IN_SUITE(TestName, SuiteName)`**,
  which expands to `TEST_F(SuiteName, TestName)`.
- Each application defines a **fast suite** fixture class
  (e.g., `KratosStructuralMechanicsFastSuite`) that inherits from
  `KratosCoreFastSuite`, creates a `Kernel`, and imports the application.
- Assertion macros live in two headers:
  - **`includes/checks.h`** — `KRATOS_CHECK_*` macros for production-code
    precondition checks (throw `KRATOS_ERROR` on failure).
  - **`includes/expect.h`** — `KRATOS_EXPECT_*` macros for test assertions
    (GTest-backed: `KRATOS_EXPECT_NEAR`, `KRATOS_EXPECT_EQ`,
    `KRATOS_EXPECT_VECTOR_NEAR`, `KRATOS_EXPECT_MATRIX_NEAR`,
    `KRATOS_EXPECT_EXCEPTION_IS_THROWN`, etc.).
- **Prefer `KRATOS_EXPECT_*` in tests** and `KRATOS_CHECK_*` in production code.
- Test files are named `test_<feature_name>.cpp` and live under
  `tests/cpp_tests/` within each application.

  ```cpp
  #include "testing/testing.h"
  #include "<app>_fast_suite.h"
  #include "containers/model.h"

  namespace Kratos::Testing {

  KRATOS_TEST_CASE_IN_SUITE(MyFeatureDoesX, KratosCoreFastSuite)
  {
      Model current_model;
      auto& r_model_part = current_model.CreateModelPart("TestPart", 1);
      // Arrange / Act
      KRATOS_EXPECT_NEAR(result, expected, 1e-10);
  }

  } // namespace Kratos::Testing
  ```

#### C++ Benchmarks (Google Benchmark)
- Benchmarks are defined using **Google Benchmark**.
- Benchmark files are named `benchmark_<feature_name>.cpp`.
- Built only when `KRATOS_BUILD_BENCHMARK=ON`.
- Always use `benchmark::State` and the `BENCHMARK()` macro.

  ```cpp
  #include <benchmark/benchmark.h>

  static void BM_MyFeature(benchmark::State& state) {
      for (auto _ : state) {
          benchmark::DoNotOptimize(/* ... */);
      }
  }

  BENCHMARK(BM_MyFeature)->Range(8, 8 << 10);
  BENCHMARK_MAIN();
  ```

---

### 🔄 GitHub Actions CI/CD

- CI workflows are defined in **`.github/workflows/`**:
  - `ci.yml` — runs on pull requests to `master`; builds and tests on
    Ubuntu (gcc + clang) with `Custom` and `FullDebug` build types.
  - `nightly_build.yml` — scheduled nightly; broader test coverage.
- When suggesting CI changes, follow the existing workflow structure.
- Do **not** suggest GitLab CI, Travis CI, or other CI systems —
  this project uses **GitHub Actions**.

---

### 🧩 Key Kratos Concepts to Be Aware Of

| Concept            | Description                                                                 |
|--------------------|-----------------------------------------------------------------------------|
| `ModelPart`        | Container for nodes, elements, conditions, and sub-model-parts              |
| `Model`            | Top-level container that owns all `ModelPart` instances                     |
| `Node`             | Geometric point with degrees of freedom (DOFs) and historical data          |
| `Element`          | Finite element — implements `CalculateLocalSystem`, etc.                    |
| `Condition`        | Boundary condition entity                                                   |
| `Process`          | Encapsulates an operation on a `ModelPart`                                  |
| `Variable`         | Typed data field (e.g., `DISPLACEMENT`, `TEMPERATURE`)                      |
| `ConstitutiveLaw`  | Material law abstraction                                                    |
| `ProcessInfo`      | Stores solver-level metadata (time step, iteration count, etc.)             |
| `DataCommunicator` | Abstraction for MPI communication                                           |
| `Kernel`           | Bootstraps the framework; loads applications                                |
| `Parameters`       | JSON-backed configuration object used for data-driven design                |

---

## Instructions

### 1) Purpose of this file

This file defines **repository-specific rules** for coding agents and
contributors in **Kratos Multiphysics**. The goal is to keep changes:

- Technically correct for the Kratos core and applications,
- Consistent with the existing Kratos code style,
- Aligned with local build/test/CI workflows,
- Minimal and safe for production and release branches.

---

### 2) Authoritative sources and precedence

When instructions conflict, use this order:

1. **Direct user request**
2. **This file (`.github/copilot-instructions.md`)**
3. **Repository code and scripts (ground truth)**
4. Generic conventions

If uncertain, prefer existing patterns in `kratos/` and `applications/`
over assumptions.

---

### 3) Project identity and boundaries

Kratos Multiphysics is the upstream framework. Both the core (`kratos/`)
and all applications (`applications/`) are first-class parts of this repository.

#### Scope guidelines

- Changes to `kratos/` (core) are valid and expected when working on
  core features, utilities, or fixes.
- Changes to `applications/` are valid for application-specific work.
- `external_libraries/` contains vendored third-party code — do **not** modify
  unless explicitly requested.
- Keep changes focused; do not refactor unrelated modules.

---

### 4) Real repository layout (important areas)

- Root build/config: `CMakeLists.txt`, `INSTALL.md`, `CONTRIBUTING.md`
- Core framework: `kratos/` (includes, sources, python, tests, benchmarks)
- Applications: `applications/*Application/`
- Configure scripts: `scripts/standard_configure.*`
- Local build scripts: `build/configure.bat`, `build/configure.sh`
- VS Code automation: `.vscode/tasks.json`
- CI workflows: `.github/workflows/`
- Runtime/install output: `bin/<BuildType>/...`
- Build trees: `build/<BuildType>/...`

#### Generated artifacts

Avoid manual edits unless explicitly requested:

- `build/<BuildType>/compile_commands.json`
- Any CMake-generated cache or install manifests

---

### 5) Build system (must follow repo workflow)

#### Primary configure/build entrypoints

- **Windows:** `scripts/standard_configure.bat` (template);
  locally, `build/configure.bat` (customized copy).
- **Linux:** `scripts/standard_configure.sh` (template);
  locally, `build/configure.sh` (customized copy).

These scripts set compilers, select applications via `KRATOS_APPLICATIONS`,
invoke CMake configure + build, and install to `bin/<BuildType>`.

#### Key build environment variables

- `KRATOS_BUILD_TYPE` (`Release`, `RelWithDebInfo`, `FullDebug`, `Custom`)
- `KRATOS_SOURCE` — path to the repository root
- `KRATOS_BUILD` — build tree root (default: `<repo>/build`)
- `KRATOS_APPLICATIONS` — semicolon-separated list of application paths
- `PYTHON_EXECUTABLE` — Python interpreter path
- `CMAKE_GENERATOR` — e.g. `Ninja`, `Visual Studio 16 2019`
- `NUMBER_OF_COMPILATION_CORES`
- `PYTHONPATH` — must include `bin/<BuildType>` at runtime
- `LD_LIBRARY_PATH` (Linux) / `PATH` (Windows) — must include `bin/<BuildType>/libs`

---

### 6) Preferred local execution path (VS Code tasks)

When running in this workspace, prefer existing tasks from `.vscode/tasks.json`:

#### Build tasks

- `Build` — configure + compile (Windows/Linux, selectable build type & generator)
- `MPI Build` — build with `KRATOS_MPI_BUILD=ON`

#### Run/test tasks

- `Run Tests` — run full Python test suite
- `Run CurrentFile` — run the currently open Python file with the correct env
- `Run C++ Tests` — run all C++ GTest suites via `run_cpp_tests.py`
- `Run C++ Test Suite` — run a specific GTest executable (e.g., `KratosCoreTest`)
- `Run C++ Test Suite Filtered` — run with `--gtest_filter=*<pattern>*`
- `Run Current Benchmark file to JSON` — run a benchmark and write JSON output

If a suitable task exists, use it instead of inventing custom command flows.

---

### 7) C++ coding conventions (repo-accurate)

Follow Kratos conventions used throughout `kratos/` and `applications/`.

#### Naming and structure

- Classes/types: `PascalCase`
- Methods (including Kratos interface overrides): `PascalCase`
  (`CalculateLocalSystem`, `GetDofList`, `Check`, etc.)
- File names: `snake_case` (e.g., `total_lagrangian.h`)
- Variables: descriptive; standard Kratos prefixes:
  - References: `rSomething`
  - Pointers: `pSomething`
  - Members: `mSomething`
- Constants and Kratos variables follow existing project style.

#### Error handling and logging

- Prefer Kratos macros in implementation code:
  - `KRATOS_TRY` / `KRATOS_CATCH("")`
  - `KRATOS_ERROR_IF`, `KRATOS_ERROR_IF_NOT`
  - `KRATOS_INFO`, `KRATOS_WARNING`
- Avoid `std::cout` for production diagnostics.

#### Memory and types

- Prefer smart pointers and Kratos pointer macros/type aliases.
- Use `const` aggressively where appropriate.

---

### 8) Python conventions

- Follow PEP 8 for new code where practical.
- Use `snake_case` for functions/variables.
- Import Kratos as `import KratosMultiphysics` and applications as
  `import KratosMultiphysics.<AppName>Application`.
- Use `KratosUnittest.TestCase` as the base for all Python tests.
- Reuse existing test harness style based on `KratosUnittest` in
  application test entry files.

---

### 9) pybind11 binding conventions

Bindings live under each application's `custom_python/` (or `kratos/python/`
for core).

- Keep binding files thin; business logic belongs in C++ classes/processes.
- Expose classes/processes with
  `py::class_<Derived, Derived::Pointer, Base>`.
- Each binding category gets its own file:
  `add_custom_<category>_to_python.cpp`.
- The main module file (`<app>_python_application.cpp`) calls all
  `AddCustom*ToPython(m)` functions and registers variables.
- Keep signatures aligned with existing modules.
- Maintain naming consistency with neighboring bindings.

Note: Prefer adding docstrings for **newly introduced public bindings**
but avoid noisy rewrites of unrelated legacy bindings.

---

### 10) CMake / application pattern to preserve

Each Kratos application follows this standard CMake pattern:

1. `file(GLOB_RECURSE ...)` — collect core sources from `custom_*/*.cpp`
2. `add_library(Kratos<App>Core SHARED ...)` — core library, links `KratosCore`
3. `pybind11_add_module(Kratos<App>Application ...)` — pybind module, links core
4. `kratos_add_gtests(TARGET Kratos<App>Core SOURCES ...)` — C++ tests
   (when `KRATOS_BUILD_TESTING=ON`)
5. Optionally add benchmarks when `KRATOS_BUILD_BENCHMARK=ON`
6. `kratos_python_install(...)` — install Python scripts/tests
7. `install(TARGETS ...)` — install shared libs into `libs/`

When adding new source files, ensure they are included by current globbing
conventions and linked correctly.

---

### 11) Testing conventions in this repo

#### C++ tests (actual pattern)

- Framework: Kratos wrappers over GTest.
- Typical includes: `testing/testing.h` and the application's
  `<app>_fast_suite.h`.
- Primary macro: `KRATOS_TEST_CASE_IN_SUITE(TestName, SuiteName)`.
- Assertion macros: `KRATOS_EXPECT_*` from `includes/expect.h`
  (e.g., `KRATOS_EXPECT_NEAR`, `KRATOS_EXPECT_EQ`, `KRATOS_EXPECT_VECTOR_NEAR`).
- Fast suite fixture classes exist per application
  (e.g., `KratosStructuralMechanicsFastSuite`) inheriting `KratosCoreFastSuite`.

#### Python tests (actual pattern)

- Application-level suite files (`test_<App>.py`) assemble
  `small`, `nightly`, `validation`, `all` using `KratosUnittest`.
- Runner scripts under `bin/<BuildType>/KratosMultiphysics/testing/`.
- Test factory classes are common for running parametric analyses.

#### MPI tests

- Use dedicated MPI test tasks/scripts.
- Keep `OMP_NUM_THREADS=1` and MPI env assumptions consistent
  with existing tasks/CI.

---

### 12) CI/CD constraints (GitHub Actions)

- CI is GitHub Actions-based (`.github/workflows/`).
- `ci.yml` — PR checks: builds Ubuntu matrix (gcc/clang × Custom/FullDebug),
  runs C++ and Python tests per changed applications.
- `nightly_build.yml` — scheduled nightly: broader test coverage including
  Windows and Rocky Linux.
- Preserve existing workflow structure and job naming conventions.
- Do **not** suggest GitLab CI, Travis CI, or other CI systems.

---

### 13) Change checklist for coding agents

Before finalizing a change:

1. Do **not** modify `external_libraries/` unless explicitly requested.
2. Keep style consistent with neighboring files.
3. Update relevant registration/binding/CMake hooks when adding new entities.
4. Run the most specific available test/task first, then broader ones if needed.
5. Avoid modifying generated/build artifacts unless requested.
6. Summarize what changed and what was validated.

---

### 14) Do / Don't summary

#### Do ✅

- Use configure scripts (`scripts/standard_configure.*` or local
  `build/configure.*`) for builds.
- Prefer existing VS Code tasks for build/test.
- Keep pybind11 files as binding glue, not logic containers.
- Make minimal, targeted changes.
- Follow Kratos naming conventions strictly in all C++ code.
- Always include header guards or `#pragma once` in C++ headers.
- Use `KRATOS_TRY` / `KRATOS_CATCH("")` macros in C++ method bodies
  for error handling where appropriate.
- Use `KRATOS_ERROR_IF` and `KRATOS_ERROR_IF_NOT` for precondition checks.
- Use `KRATOS_INFO`, `KRATOS_WARNING` for logging — never use `std::cout`
  directly in production code.
- Use `KRATOS_EXPECT_*` macros in C++ tests; use `KRATOS_CHECK_*` in
  production-code preconditions.
- Write docstrings/comments for all public APIs (both C++ and Python).
- Keep pybind11 binding files thin — logic belongs in C++, not in binding code.
- Prefer `const auto&` for range-based loops over containers.
- Always initialize variables at declaration.

#### Don't ❌

- Do **not** add unrelated refactors while fixing targeted issues.
- Do **not** modify files inside `external_libraries/` by default.
- Do **not** use raw `new`/`delete` — use smart pointers.
- Do **not** use `std::cout` for logging in production C++ code.
- Do **not** hardcode file paths — use Kratos path utilities or CMake variables.
- Do **not** use `using namespace std;` globally.
- Do **not** suggest GitLab CI, Travis CI, or other CI systems —
  this project uses **GitHub Actions**.
- Do **not** bypass the configure scripts with ad-hoc CMake commands.
- Do **not** write Python tests without inheriting from
  `KratosUnittest.TestCase`.
- Do **not** hand-edit generated build artifacts.

---

### 15) Additional notes

- When generating new **Kratos applications**, follow the standard application
  scaffolding structure exactly (see `applications/StructuralMechanicsApplication/`
  as a reference).
- When generating new **Elements or Conditions**, always override the required
  virtual methods: `CalculateLocalSystem`, `EquationIdVector`,
  `GetDofList`, and `Check`.
- When generating new **Processes**, inherit from `Process` and implement
  `Execute()` or the appropriate stage methods (`ExecuteInitialize`,
  `ExecuteBeforeSolutionLoop`, `ExecuteInitializeSolutionStep`, etc.).
- Prefer **data-driven design** using Kratos `Parameters` (JSON-based) for
  process and solver configuration.
- When writing pybind11 bindings for Processes, expose the `Execute()` method
  and any stage-specific methods.
- Branch naming convention: `subject/short-description`
  (e.g., `core/adding-xxx-utility`, `core/fix-xxx-utility`).

---

### 16) Notes for future instruction maintenance

If repository workflows change, update this file using evidence from:

- root `CMakeLists.txt` and `INSTALL.md`,
- configure scripts in `scripts/`,
- `.vscode/tasks.json`,
- `.github/workflows/*.yml`,
- representative files under `kratos/` and `applications/`.

This file should remain practical and repo-grounded, not generic.
