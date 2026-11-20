---
title: Testing and contributing
keywords: tests suite small nightly validation mpi ci torch-free guards contributing checklist mutation
tags: [Testing_And_Contributing.md]
sidebar: physicsnemo_application
summary: How the test suite is organised and run, what the torch-free CI actually exercises, the five guards that keep the suite and the documentation honest, and the checklist for adding a module, a process or a page.
---

# Testing and contributing

## The suites

| Suite | What it holds | Where it runs |
|---|---|---|
| `small` | everything serial - 882 tests here, under a minute; ML tests self-skip without their package | the GitHub CI (torch-free), every local run |
| `nightly` | the small suite | the nightly workflow |
| `validation` | the small suite plus the 19 example notebooks executed end to end (about five minutes) | on request only |
| the MPI runner `test_PhysicsNeMoApplication_mpi.py` | 37 tests over gloo at two and three ranks, no Metis required | on request |

```bash
cd applications/PhysicsNeMoApplication/tests
python3 test_PhysicsNeMoApplication.py -l small            # or nightly, validation
python3 test_graph_bridge.py                               # one file
python3 test_graph_bridge.py TestBuildGraph.test_EdgesAreBidirectional   # one test
OMP_NUM_THREADS=1 mpiexec -np 3 python3 test_PhysicsNeMoApplication_mpi.py --using-mpi
```

`tests/` mirrors `python_scripts/` by name (`test_<module>.py`), with subdirectories for `bridges/mesh_bridge/` and `active_learning/`. The real-solver cases live in `tests/kratos_solver_cases/` (thermal, structural, fluid, thermo-mechanical, a transient harness) and `tests/adjoint_cases/`; they are built in memory, without `.mdpa` files, and are what the notebooks and the Examples cases import too.

**Two development gotchas.** `bin/Release/KratosMultiphysics/PhysicsNeMoApplication/` holds *copies* of `python_scripts/` made at install time - after editing a module, copy it there (or rebuild) before running tests, or a stale copy silently runs the old code. And small solves are much faster with `OMP_NUM_THREADS=2` than with the default twenty threads.

## What the torch-free CI exercises

The Linux CI builds the application and runs the small suite on a runner with no torch, no physicsnemo and no GPU. Every ML-dependent test class is gated on a lazy import and skips there; what remains - and it is a lot - is the export processes, the tessellation and provenance machinery, the graph and grid bridges' numpy paths, the settings validation of every process, the import contract, the documentation guards and the benchmark smoke test. Locally you can reproduce that run by putting a directory on `PYTHONPATH` that holds shim packages (`torch/__init__.py`, `physicsnemo/__init__.py`, ...) whose only content is `raise ImportError`; the suite then reports about 450 skips and must still pass.

The lesson that made this policy: when upgrading physicsnemo, **compare skip counts, not pass/fail**. A renamed upstream symbol turns whole test classes into green skips.

## The five guards

`tests/test_suite_registration.py` is torch-free, `ast`-only, and fails the suite when any of these rot:

1. **Registration.** Every `TestCase` subclass under `tests/` (subdirectories included) must appear inside a `loadTestsFromTestCases([...])` call in a runner. Importing a class is not registering it; an unregistered class runs zero tests and the count still looks healthy.
2. **Documented identifiers.** Every backticked `module.Attribute` in the README and the documentation pages, where `module` is an application module, must name something that module defines.
3. **Documented import paths.** Every written `KratosMultiphysics.PhysicsNeMoApplication.<path>` must resolve in the tree, and a trailing attribute must be defined by the module it names.
4. **Documented process factories.** Every `"python_module"`/`"kratos_module"` pair in a documented JSON block must resolve to a module with a `Factory`.
5. **Cross-module attributes.** `from ... import mod` followed by `mod.Name` in any source or test must reach a name `mod` still defines - the failure that splitting a module produces without breaking a single import.

Plus `tests/test_import_contract.py` (the application imports with no ML package installed), `tests/test_docs_links.py` (every relative link and image reference in the documentation resolves) and `tests/test_docs_figures.py` (the figure generator still produces its files). The benchmark script is executed by a smoke test so its quoted numbers cannot drift unnoticed.

What the guards deliberately cannot catch: **stale prose about project state**. A sentence saying a feature is "not yet shipped" passes every identifier check after the feature ships. Status prose has to be read against the tree, which is what the roadmap rounds do.

## Writing tests that would fail on a wrong answer

Several real bugs shipped through a green suite before these rules were written down:

- **Stand-in models must be affine, not linear.** A model with `f(0) = 0` cannot detect a missing de-normalization *shift*; three de-normalization bugs survived because every stand-in was a `Doubler`.
- **Fixtures must not be a unit cube.** On `[0, 1]^3` every length scale is one and `normalize_coordinates` is the identity; the shared fixture takes an `extent` now.
- **Canaries must be asymmetric.** A gradient-layout flip upstream passed 469 tests because the test field's gradient was a symmetric matrix.
- **Single-step tests cannot see per-step waste.** A process that rebuilt the whole mesh graph every step (1073 ms to keep 2.8 ms of features) was invisible until a test ran it twice.
- **Never validate a refactored helper against the function that now calls it** - both sides move together. Validate against the mesh.
- **Collectives must run on all ranks.** A collective inside a rank-0 guard deadlocks - sometimes only after an unrelated import changes timing.
- **Mutation-check every pin.** Delete the call the test protects; the test must fail. Then confirm the same test id passes on clean code, because `FAILED (errors=1)` from a wrong test name is not a kill.

## Adding things

**A module.** Pick the package by the rule in `python_scripts/README.md` (a `Factory` means `processes/`; conversion means `bridges/`; and so on), keep every `import torch` / `import physicsnemo` / `import cupy` inside a lazy helper, add `tests/test_<module>.py`, register its classes in the runner, and name the module on [Where things live](Module_Map.html).

**A process.** Subclass `Kratos.Process` (or `InferenceProcess` for anything that runs a model - the card, normalization, guard and uncertainty machinery come for free), give it a `Factory(settings, model)`, validate settings with a default block, extract topology in `ExecuteInitialize` and only values per step, and add it to the [Process reference](Process_Reference.html) with a JSON block - the guard resolves it.

**A documentation page.** Front matter with `title`, `keywords`, `tags`, `sidebar: physicsnemo_application` and a one-line `summary` with **no colons** in any value (the site generator splits on them); a category needs a `menu_info.json` listing its pages in order and an entry in the application's top-level `menu_info.json`; images go in the category's `images/` folder - a diagram as an SVG source with its PNG twin next to it (`make_docs_figures.py --rasterize` writes the twin; both are committed), a data-driven figure as a PNG from the generator; mermaid diagrams use the `<div class="mermaid">` form the site renders; only backtick `module.Attribute` when both exist.

**A notebook.** Number it, add a row to `examples/notebooks/README.md` and to [Examples](../Examples/Examples.html); `tests/test_notebooks.py` picks it up automatically and runs it in a throwaway copy of the tree on the validation suite.

## Branch and commit conventions

The Kratos branch convention is `subject/short-description`; this application's history uses a `[PhysicsNeMo]` prefix on commit subjects. Do not modify `external_libraries/`, do not add the CMake gate, and keep the CI JSON lists byte-exact when touching them (they have no trailing newline).

Next: [Troubleshooting and traps](Troubleshooting_And_Traps.html).

## Solver cases, and the fixtures they borrow

`tests/kratos_solver_cases/` holds the real solves the bridges are tested against. The older ones are built entirely IN MEMORY - the thermal, structural, fluid-cavity and coupled sintering cases construct their mesh and add every solution-step variable before the nodes exist, because `use_input_model_part` fails at construction otherwise.

The cases added for the compiled-but-unexercised applications work differently. Each one REUSES the owning application's own test data in place rather than copying it, and reports its absence instead of failing:

| Case | Drives | Cost on the reference machine |
|---|---|---|
| `cylinder_case` | transient flow past a cylinder (FluidDynamics) | 1110 nodes, 0.2 s to t = 0.2 |
| `rans_case` | k-epsilon channel flow (RANS) | 49 nodes, 0.04 s |
| `airfoil_case` | NACA 0012 in compressible potential flow | 22 nodes, 0.02 s |
| `geomechanics_case` | one-dimensional consolidation (GeoMechanics) | 201 nodes, 0.14 s per stage |
| `contact_case` | an ALM frictionless patch test (ContactStructuralMechanics) | 8 nodes, 0.05 s |
| `iga_case` | the Scordelis-Lo roof (Iga) | 36 control points |
| `compliance_case` | a SIMP cantilever, built in memory | scales with the design grid |

Borrowing rather than copying keeps a 216 kB cylinder mesh out of this application and means the fixtures cannot drift from the cases they came from. Each module exposes `IsAvailable()`, and its tests skip when the owning application's data is not there.

What borrowing costs is a little parameter surgery, and the traps are worth knowing:

- **Reference-file comparisons travel with the parameters.** Every borrowed case carries processes that compare against `.json` or `.dat` files relative to the owning directory. They belong to that application's regression, so each case filters them out BY PROCESS NAME - the lists they sit in are named differently from one application to the next (`json_check_process` here, `json_check_process_list` there).
- **Dropping the whole `processes` block is not the shortcut it looks like.** It also removes the boundary conditions: the cylinder then runs to completion in no time at all with a maximum velocity of exactly zero, which reads like a fast case rather than an undriven one.
- **Paths hide in more places than `model_import_settings`.** The contact mesh is named inside a MODELER's own import settings, two levels down, and as an ARRAY of two blocks. The materials file, the geometry and physics files of an IGA modeler, and the results files all need the same treatment.
- **The RANS parameters are a template**, carrying placeholders like `<STABILIZATION_METHOD>` that the owning driver substitutes per test; the case substitutes them rather than shipping a second copy.
- **A case that runs is not a case that solves.** IgaApplication's `single_patch_test` is structurally singular - a direct solver reports a zero column, and the iterative solver it is configured with returns NaN displacements without complaining - which its own test never notices because it asserts nothing numerically. The Scordelis roof is used instead.
