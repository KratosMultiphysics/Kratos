# AGENTS.md — Kratos Multiphysics

Entry point for any coding agent working in this repository. Read this file first, then the
path-scoped conventions indexed at the bottom.

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
├── AGENTS.md                          # This file — project overview and global rules
├── .github/
│   ├── copilot-instructions.md        # Thin pointer to /AGENTS.md
│   ├── instructions/                  # Path-scoped convention files (see index below)
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
├── docs/                              # Documentation sources
├── documents/                         # Additional project documents
├── CMakeLists.txt                     # Root CMake entry point
├── README.md                          # Project introduction
├── INSTALL.md                         # Build instructions
├── CONTRIBUTING.md                    # Contribution guidelines
├── CODE_OF_CONDUCT.md                 # Community guidelines
└── license.txt                        # BSD license
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
2. **`AGENTS.md` and `.github/instructions/` files**
3. **Repository code and scripts (ground truth)**
4. Generic conventions

If uncertain, prefer existing patterns in `kratos/` and `applications/` over assumptions.

---

## Project Boundaries

- Changes to `kratos/` (core) and `applications/` are both valid and expected.
- `external_libraries/` contains vendored third-party code — do **not** modify unless explicitly requested.
- Keep changes focused; do not refactor unrelated modules.

---

## Agent Behavior Principles

### 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

- **State assumptions explicitly** — if uncertain about intent, ask rather than guess.
- **Present multiple interpretations** — when ambiguity exists, name the options; don't pick silently.
- **Push back when warranted** — if a simpler approach exists, say so before implementing.
- **Stop when confused** — name what's unclear and ask for clarification rather than running with a wrong assumption.

### 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was explicitly asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for scenarios that cannot happen.
- If 200 lines could be 50, prefer the 50-line version.

**The test:** Would a senior engineer say this is overcomplicated? If yes, simplify.

### 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting that wasn't part of the request.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it — don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that **your** changes made unused.
- Don't remove pre-existing dead code unless explicitly asked.

**The test:** Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform imperative tasks into verifiable goals:

| Instead of... | Transform to... |
|---------------|-----------------|
| "Add validation" | "Write tests for invalid inputs, then make them pass" |
| "Fix the bug" | "Write a test that reproduces it, then make it pass" |
| "Refactor X" | "Ensure tests pass before and after" |

For multi-step tasks, state a brief plan before starting:

```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

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

Topic-specific rules live in `.github/instructions/`. **Read the matching file before editing
files of that kind.** GitHub Copilot loads these automatically from their `applyTo` globs; other
agents should open them on demand using this index.

| File | Read when working on |
|------|----------------------|
| `build.instructions.md` | Building, running tests, configure scripts, VS Code tasks |
| `ci-cd.instructions.md` | `.github/workflows/*.yml` |
| `cmake.instructions.md` | `**/CMakeLists.txt` |
| `cpp-conventions.instructions.md` | `**/*.cpp`, `**/*.h`, `**/*.hpp` |
| `python-conventions.instructions.md` | `**/*.py`, `**/custom_python/**/*.cpp` |
| `testing.instructions.md` | `**/tests/**`, `test_*.cpp`, `test_*.py`, `benchmark_*.cpp` |

---

## Tool-Specific Extras

Optional helpers, written in the GitHub Copilot file formats:

| Path | Contents |
|------|----------|
| `.github/agents/` | Custom agent personas (e.g. `kratos-reviewer` — read-only convention reviewer) |
| `.github/prompts/` | Reusable task prompts (scaffold a `Process`, scaffold a test) |
| `.github/skills/` | Multi-step workflows (e.g. `scaffold-application`) |

Agents that don't understand these formats can still read the files as plain Markdown and ignore
the YAML frontmatter — the instructions inside are tool-neutral.
