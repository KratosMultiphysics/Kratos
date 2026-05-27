# GitHub Copilot Instructions — Kratos Multiphysics

## Project Overview

**Kratos Multiphysics** is an open-source framework for building parallel,
multi-disciplinary simulation software. It provides a core finite-element
engine (`kratos/`) plus a library of domain-specific **applications**
(structural mechanics, fluid dynamics, contact, etc.) and exposes its
full API to Python via **pybind11**.

---

## Repository Structure

```
.
├── .github/
│   ├── copilot-instructions.md        # This file — project overview and global rules
│   ├── instructions/                  # Scoped instruction files (auto-loaded per file type)
│   ├── prompts/                       # Reusable task prompts
│   ├── agents/                        # Custom agent personas
│   ├── skills/                        # On-demand multi-step workflow skills
│   └── workflows/                     # GitHub Actions CI definitions
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
│       ├── <application_name>_application.h/.cpp
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
├── external_libraries/                # Vendored third-party deps (amgcl, …)
├── scripts/                           # Configure script templates
│   ├── standard_configure.sh
│   ├── standard_configure.bat
│   └── ...
├── cmake_modules/                     # Custom CMake Find/utility modules
├── CMakeLists.txt                     # Root CMake entry point
├── INSTALL.md                         # Build instructions
└── CONTRIBUTING.md                    # Contribution guidelines
```

---

## Tech Stack

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
| Dependencies     | Bundled in `external_libraries/`                             |

---

## Key Kratos Concepts

| Concept            | Description                                                                 |
|--------------------|-----------------------------------------------------------------------------|
| `Model`            | Top-level container that owns all `ModelPart` instances                     |
| `ModelPart`        | Container for nodes, elements, conditions, and sub-model-parts              |
| `Node`             | Geometric point with degrees of freedom (DOFs) and historical data          |
| `Element`          | Finite element — implements `CalculateLocalSystem`, etc.                    |
| `Condition`        | Boundary condition entity                                                   |
| `Process`          | Encapsulates an operation on a `ModelPart`                                  |
| `Variable`         | Typed data field (e.g., `DISPLACEMENT`, `TEMPERATURE`)                      |
| `ConstitutiveLaw`  | Material law abstraction                                                    |
| `ProcessInfo`      | Stores solver-level metadata (time step, iteration count, etc.)             |
| `Parameters`       | JSON-backed configuration object used for data-driven design                |
| `Kernel`           | Bootstraps the framework; loads applications                                |
| `DataCommunicator` | Abstraction for MPI communication                                           |
| `Strategy`         | Top-level solver orchestration                                              |
| `Scheme`           | Time integration / linearization                                            |
| `BuilderAndSolver` | Assembles the global system and solves                                      |

---

## Architecture

- **`kratos/`** is the framework core: containers, solvers, geometries, I/O, and the Python binding layer.
- **`applications/`** extend the core with domain-specific Elements, Conditions, ConstitutiveLaws, Processes, and Strategies — without modifying the core.
- Each application registers its components via `RegisterComponents()` in `<app>_application.cpp`.
- Branch naming convention: `subject/short-description` (e.g., `core/adding-xxx-utility`, `structural/fix-xxx-element`).

---

## Authoritative Sources and Precedence

When instructions conflict, use this order:

1. **Direct user request**
2. **This file and `.github/instructions/` files**
3. **Repository code and scripts (ground truth)**
4. Generic conventions

If uncertain, prefer existing patterns in `kratos/` and `applications/` over assumptions.

---

## Project Boundaries

- Changes to `kratos/` (core) and `applications/` are both valid and expected.
- `external_libraries/` contains vendored third-party code — do **not** modify unless explicitly requested.
- Keep changes focused; do not refactor unrelated modules.

---

## Change Checklist

Before finalizing a change:

1. Do **not** modify `external_libraries/` unless explicitly requested.
2. Keep style consistent with neighboring files.
3. Update CMake, variable registration, and pybind11 binding hooks when adding new entities.
4. Register new components in `RegisterComponents()` inside `<app>_application.cpp`.
5. Run the most specific available test/task first, then broader ones if needed.
6. Avoid modifying generated/build artifacts unless requested.

---

## Detailed Conventions

Topic-specific rules live in `.github/instructions/`:

| File | Applies to |
|------|-----------|
| `agent-behavior.instructions.md` | All files — always active |
| `build.instructions.md` | Build scripts, configure scripts, VS Code tasks |
| `ci-cd.instructions.md` | `.github/workflows/*.yml` |
| `cmake.instructions.md` | `CMakeLists.txt` files |
| `cpp-conventions.instructions.md` | `*.cpp`, `*.h`, `*.hpp` files |
| `python-conventions.instructions.md` | `*.py` and `custom_python/**/*.cpp` binding files |
| `testing.instructions.md` | `tests/` directories, test and benchmark files |
