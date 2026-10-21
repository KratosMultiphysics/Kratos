# AGENTS.md: GeoMechanicsApplication

Entry point for any coding agent that works in `applications/GeoMechanicsApplication/`.
Read the root `/AGENTS.md` first, then this file. Inside this folder, the rules here and in
`agents/` (next to this file) take precedence over the general rules in the root `agents/`
folder. Where the two disagree, this file says so explicitly.

Maintainers: `@KratosMultiphysics/geomechanics` (see `.github/CODEOWNERS`). Questions about
undocumented behaviour go to them.

## What the application does

Geotechnical analyses: staged analysis, K0 procedure, quasi-static and dynamic U-Pw
(displacement + water pressure) analysis, groundwater flow (Pw), heat transport (T),
interface elements, trusses, cables, beams, backward erosion piping, and user-defined soil
models loaded from UDSM (PLAXIS format) and UMAT (ABAQUS format) libraries. See `README.md`.

The application depends on `StructuralMechanicsApplication` and `LinearSolversApplication`.
MPI is not supported.

## Folder map

| Path | Contents |
|------|----------|
| `custom_elements/` | Elements. `contribution_calculators/` holds reusable LHS/RHS calculators; stress state and integration coefficient policies are described in `Policies.md` |
| `custom_conditions/` | Boundary conditions (loads, fluxes, Lysmer absorbing, micro-climate, ...) |
| `custom_constitutive/` | Constitutive laws, UDSM/UMAT wrappers, Mohr-Coulomb, cap models |
| `custom_retention/` | Retention laws (`SaturatedLaw`, `VanGenuchtenLaw`, ...) |
| `custom_geometries/` | Line and surface interface geometries (`interface_geometry.hpp`) |
| `custom_processes/` | C++ processes |
| `custom_strategies/` | Schemes, strategies, builders and solvers |
| `custom_utilities/` | Utility classes (static member functions) and factories |
| `custom_workflows/` | C++ analysis routes (`dgeosettlement`, `dgeoflow`) with their own process factory |
| `custom_operations/` | Kratos operations |
| `custom_python/` | pybind11 bindings |
| `python_scripts/` | Analysis stage, solvers, Python process wrappers |
| `test_setup_utilities/` | Shared C++ helpers for unit tests and benchmarks (`ModelSetupUtilities`, `ElementSetupUtilities`) |
| `tests/cpp_tests/` | C++ unit tests; sub-folders mirror the `custom_*` folders |
| `tests/` | Python integration, regression and validation tests with their input data |
| `benchmarks/` | Google Benchmark executables |
| `documents/` | LaTeX verification report and tutorials |

Most `custom_*` folders have a `README.md` with the theory and the design of their content.
Read the README of the folder you change, and update it when you change documented behaviour.

## Agent files for this application

| File | Read when |
|------|-----------|
| `agents/instructions/geo-cpp.instructions.md` | Editing any C++ file in this application |
| `agents/instructions/geo-testing.instructions.md` | Writing or running C++ or Python tests |
| `agents/instructions/geo-python.instructions.md` | Editing Python files or pybind11 bindings |
| `agents/instructions/geo-domain.instructions.md` | Touching elements, conditions, laws, processes or their physics |
| `agents/reviewers/geo-reviewer.agent.md` | Reviewing a change in this application |
| `agents/skills/geo-build-and-test/SKILL.md` | Building the application and running its tests |
| `agents/skills/geo-add-process/SKILL.md` | Adding a new process |
| `agents/prompts/geo-scaffold-test.prompt.md` | Creating a new C++ or Python test |

## Change checklist

1. Run `clang-format` with `applications/GeoMechanicsApplication/.clang-format` on every C++
   file you change.
2. Add a unit test in `tests/cpp_tests/<same sub-folder as the source>/` for new C++ code.
3. New source files in an existing `custom_*` folder are picked up by `GLOB_RECURSE` in
   `CMakeLists.txt`; a new top-level folder must be added to that glob.
4. Register new elements, conditions and constitutive laws in
   `geo_mechanics_application.h/.cpp`, and new variables in
   `geo_mechanics_application_variables.h/.cpp` and `custom_python/geo_mechanics_python_application.cpp`.
5. Add new Python test classes to the right suite in `tests/test_GeoMechanicsApplication.py`.
6. Do not add new JSON files to `exceptions_for_json_formatting.txt` (see `geo-testing.instructions.md`).


