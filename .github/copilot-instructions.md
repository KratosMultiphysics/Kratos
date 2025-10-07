# Copilot instructions for this Kratos workspace    #####CKECKLEO####

Use this to be productive quickly in this Kratos Multiphysics repo with local build/test wiring.

## Big picture
- C++ core + Python interface; applications are plugins in `applications/*` built selectively.
- Python “analyses” derive from `kratos/python_scripts/analysis_stage.py::AnalysisStage` and are configured by JSON ProjectParameters.
- Factories: `process_factory` builds Processes; `model_parameters_factory` builds Modelers from Parameters.

## Lifecycle and config
- Analysis flow: Initialize → RunSolutionLoop (InitializeSolutionStep → Predict → SolveSolutionStep → FinalizeSolutionStep → OutputSolutionStep) → Finalize.
- Modelers (from `project_parameters["modelers"]`) run before solver import: `SetupGeometryModel` → `PrepareGeometryModel` → `SetupModelPart`.
- Processes are read from `project_parameters["processes"]`; output processes must be under `output_processes` (deprecated otherwise) — see `analysis_stage.py`.
- `problem_data` must include `{echo_level, parallel_type, start_time, end_time}`; base class warns if `parallel_type` ≠ runtime (OpenMP vs MPI).

## Build (custom local script)
- Use VS Code task “Build” → `scripts/configure.sh`:
  - Activates venv named `kratos_<compiler>_<BuildType>` (e.g., `kratos_gcc_Debug`).
  - Sets `KRATOS_CPP_CONFIG_NAME=<compiler>_<BuildType>`, builds in `build/${KRATOS_CPP_CONFIG_NAME}`.
  - Installs to active site-packages (`-DCMAKE_INSTALL_PREFIX`), enables `-DKRATOS_GENERATE_PYTHON_STUBS=ON`.
  - Compiles selected apps: Iga, StructuralMechanics, Optimization, LinearSolvers, ConvectionDiffusion, ShapeOptimization, SystemIdentification.
  - Symlinks `compile_commands.json`; links installed `test/` to `build/test`.
- To change apps/flags, edit `scripts/configure.sh` (adjust `add_app` and CMake args).

## Run and test
- Task “Run all tests”: `kratos/python_scripts/testing/run_tests.py` (env `OMP_NUM_THREADS=30`). Flags: `-l small|nightly|all`, `-v 0|1|2`, `-a 'App1:App2'`.
- Task “Run cpp pattern tests”: filters C++ gtests via `.vscode/cpp_test.py` (depends on Build).
- Task “Run current file”: runs the active Python file in its folder.

## Conventions
- Analyses import `AnalysisStage`, delegate to `_GetSolver()`, and override only needed hooks.
- Processes call order: `ExecuteInitialize` → `ExecuteBeforeSolutionLoop` → per-step init/finalize → `ExecuteBeforeOutputStep`/`ExecuteAfterOutputStep` → `ExecuteFinalize`.
- Output processes go in `output_processes` and are appended to the regular list by the base stage.

## Pointers
- Base orchestration: `kratos/python_scripts/analysis_stage.py`.
- Examples: `applications/**/python_scripts/*_analysis.py` (e.g., Structural/Fluid/GeoMechanics).
- Test docs: `docs/pages/Kratos/For_Developers/CICD/Unitary_Tests.md`.

Tips: If you add a new app, wire it in `scripts/configure.sh` via `add_app`. Keep stage/solver sequencing intact. Match `parallel_type` to runtime.