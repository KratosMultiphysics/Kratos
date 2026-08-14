# Bofang initialization regression

Investigation of whether the June 2025 change to
`DamBofangConditionTemperatureProcess` (PR #13472) can modify the results or
initialization state of a thermomechanical DamApplication simulation.

The case is small, deterministic and reproducible. It establishes a clean
baseline before the DamApplication -> StructuralMechanicsApplication element
migration (handled in a separate branch).

## Historical reference

PR #13472 ("Dam/assign scalar variable process", merged 2025-06-02) changed
several Dam processes to the generic Core scalar-assignment lifecycle. For
Bofang the changes were:

| aspect            | before (ece5cfe...)              | after (359d1ef...)              |
|-------------------|----------------------------------|----------------------------------|
| fixity parameter  | `is_fixed`                       | `constrained`                    |
| lifecycle method  | `ExecuteInitialize()`            | `ExecuteBeforeSolutionLoop()`    |

Merge commit: `359d1ef414fabf1a917933cae8222b9dc657b92f`
Immediately-preceding parent: `ece5cfed146f7fcde48681c6759b3af29aa81fb0`

Confirmed with:

```
git diff ece5cfed146f7fcde48681c6759b3af29aa81fb0 359d1ef414fabf1a917933cae8222b9dc657b92f \
    -- applications/DamApplication/custom_processes/dam_bofang_condition_temperature_process.hpp
```

which shows exactly three hunks: `is_fixed` -> `constrained` (defaults and
member read) and `ExecuteInitialize()` -> `ExecuteBeforeSolutionLoop()`.

## Variants

| variant | code                                                              | purpose                      |
|---------|-------------------------------------------------------------------|------------------------------|
| A       | current master (`dam/initialization-regression` HEAD)             | reference                    |
| B       | current master + Bofang process reverted to `ExecuteInitialize()` | legacy lifecycle (process-only) |
| B2      | current master + Bofang reverted to `ExecuteInitialize()` **and** the `dam_analysis.py` pre-solver-initialize process call reverted to `ExecuteInitialize()` | faithful pre-#13472 lifecycle |
| C       | historical Kratos at `ece5cfed146` (worktree)                     | real pre-PR code             |

### Why both B and B2?

`DamAnalysis::Initialize()` is driven by `dam_analysis.py`. Before #13472 it
called `process.ExecuteInitialize()` on the processes *before*
`solver.Initialize()`; after #13472 it calls `process.ExecuteBeforeSolutionLoop()`.
The PR changed **both** the analysis call site and the Bofang process. Variant B
reverts only the process (literal reading of the task), which makes Bofang's
`ExecuteInitialize()` dead code under the new analysis (temperature is then only
assigned at `ExecuteInitializeSolutionStep`). Variant B2 reverts process **and**
analysis call site, reproducing the exact pre-PR lifecycle.

The A vs B2 comparison is the faithful isolation of the lifecycle change; A vs B
shows what happens if only the process is reverted.

## Case

2D plane-strain rectangular dam section (10 m x 20 m), 18 nodes, 10 quads.

* mechanical: `SmallDisplacementThermoMechanicElement2D4N` +
  `ThermalLinearElastic2DPlaneStrain`
* thermal: `EulerianConvDiff2D4N`
* upstream face model part `BOFANGTEMPERATURE_Reservoir_Auto1` (X = 0)
* water level 14 m covers nodes at Y = 0, 4, 8, 12 m
* Bofang parameters (fixed scalars, no tables):
  `Surface_Temp=10`, `Bottom_Temp=5`, `Height_Dam=20`, `Temperature_Amplitude=3`,
  `Day_Max_Temp=15`, `Water_level=14`, `Month=6.5`, `Gravity_Direction="Y"`,
  `Reservoir_Bottom_Coordinate_in_Gravity_Direction=0`, `constrained=false`
* materials: `DENSITY=2400`, `E=26 GPa`, `nu=0.2`, `alpha=1e-5`,
  `k=2.2`, `cp=1000`
* mechanical BC: base (Y=0) fully fixed
* solver: `dam_thermo_mechanic_solver`, direct deterministic
  `skyline_lu_factorization`, 1 thread, `theta=1.0`, `dt=86400 s`,
  `end=259200 s` (3 steps), `buffer_size=2`

The Bofang temperature varies with depth (about 4.69 C at the base to 7.29 C at
the water surface), so the initial condition is genuinely non-uniform.

> Note: `bofang_factory.py` instantiates the raw C++
> `DamBofangConditionTemperatureProcess` directly instead of the production
> `impose_reservoir_temperature_condition_process` wrapper, because the wrapper
> does not forward `ExecuteInitialize()` and would mask the lifecycle behaviour
> under study.

## Layout

```
case/                  .mdpa, ProjectParameters.json (master), bofang_factory.py
instrumentation/       instrumented_dam_analysis.py, bofang_lifecycle_recorder.py
comparison/            compare_results.py, analytical_bofang.py
results/
  current_master/                  variant A
  current_master_legacy_bofang/    variant B
  current_master_legacy_bofang_faithful/  variant B2
  pre_13472/                       variant C
run_experiment.sh      variants A, B, B2 (master build)
```

## Running

Variant A / B / B2 (from the master tree):

```
./run_experiment.sh A  <kratos_bin> <build_dir> <patch_dir>
./run_experiment.sh B  <kratos_bin> <build_dir> <patch_dir>
./run_experiment.sh B2 <kratos_bin> <build_dir> <patch_dir>
```

The patches live in `scripts/` too (see `bofang_legacy_process_only.patch`,
`bofang_legacy_faithful.patch`).

Variant C (from the historical worktree at `ece5cfed146`):

```
./run_experiment_c.sh
```

### Build commands

All builds used at most 8 parallel jobs (`-j8`).

```
# master (variants A/B/B2)
export KRATOS_APPLICATIONS="applications/StructuralMechanicsApplication;applications/PoromechanicsApplication;applications/ConvectionDiffusionApplication;applications/DamApplication"
cmake -H. -Bbuild/master -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=/usr/bin/python3 \
      -DUSE_MPI=OFF -DUSE_TRIANGLE_NONFREE_TPL=OFF -DKRATOS_BUILD_TESTING=OFF
cmake --build build/master --target install -- -j8

# variant B/B2 patch rebuild (only the DamApplication pybind module)
cmake --build build/master --target KratosDamApplication install -- -j8

# historical worktree
cd worktrees/pre_13472
export KRATOS_APPLICATIONS=...
cmake -H. -Bbuild/pre_13472 -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=/usr/bin/python3 \
      -DUSE_MPI=OFF -DUSE_TRIANGLE_NONFREE_TPL=OFF -DKRATOS_BUILD_TESTING=OFF
cmake --build build/pre_13472 --target install -- -j8
```

## Comparison

```
cd comparison
python3 compare_results.py --base-dir ../results
python3 analytical_bofang.py
```

## Findings

Experiment executed on branch `dam/initialization-regression` at
`f4b2e970ede` (master-derived HEAD), with variants A, B, B2 and C.

### Lifecycle observation (instrumented run)

The instrumented analysis records nodal `TEMPERATURE` at each lifecycle stage
(`lifecycle_temperatures.csv`). For the current implementation (A) the Bofang
value is assigned before `solver.Initialize()`:

| stage | node 7 (y=12, depth 2 m) |
|-------|--------------------------|
| after process construction | 0.0 |
| after pre-solver-init process call | 6.839771675798264 |
| before solver.Initialize | 6.839771675798264 |
| after solver.Initialize | 6.839771675798264 |
| after post-solver-init ExecuteBeforeSolutionLoop | 6.839771675798264 |

`solver.Initialize()` neither consumes nor modifies `TEMPERATURE`.

### Source-code analysis

* `SolidElement::Initialize()` only calls `InitializeConstitutiveLaw()`
  (`solid_element.cpp:352`); the constitutive law `InitializeMaterial` does not
  read `TEMPERATURE`.
* `ThermalLinearElastic3DLaw::CalculateMaterialResponseKirchhoff` is the only
  place that reads `TEMPERATURE` / `NODAL_REFERENCE_TEMPERATURE`, and it runs
  during the mechanical element's `CalculateLocalSystem`, i.e. during the solve,
  after the Bofang `ExecuteInitializeSolutionStep` of the current step.
* The Dam `IncrementalUpdateStaticSmoothingScheme::Initialize` only clears
  `INITIAL_STRESS_TENSOR`.

Consequence: whether Bofang assigns temperature at `ExecuteInitialize()` or at
`ExecuteBeforeSolutionLoop()` is numerically irrelevant as long as the value is
present before the first `solver.Solve()`.

### Numerical results (direct solver `skyline_lu_factorization`, 1 thread)

Max absolute differences:

| comparison | temperature [K] | disp. x [m] | disp. y [m] | stress xx [Pa] |
|------------|-----------------|-------------|-------------|----------------|
| A vs B2    | 0.0             | 0.0         | 0.0         | 0.0            |
| A vs C     | 0.0             | 0.0         | 0.0         | 0.0            |
| B2 vs C    | 0.0             | 0.0         | 0.0         | 0.0            |
| A vs B     | 6.84            | 1.26e-3     | 7.81e-4     | 1.68e6         |

* **A vs B2** (the faithful isolation of the PR #13472 lifecycle change) is
  bit-identical. The June 2025 change does **not** modify the results or the
  initialization state of the thermomechanical simulation.
* **A vs C** (current master vs the actual pre-PR Kratos at `ece5cfe`) is
  bit-identical: the Dam elements/constitutive laws/strategies are unchanged
  between the two commits, so the whole thermomechanical response is identical.
* **A vs B** (reverting only the C++ process to `ExecuteInitialize()` while
  keeping the post-PR `dam_analysis.py`) *does* differ: the current analysis
  never calls `ExecuteInitialize()` on its processes, so Bofang's assignment
  becomes dead code during initialization, the initial temperature is 0, and the
  whole thermomechanical response collapses to ~0. This documents the
  lifecycle coupling between `dam_analysis.py` and the Bofang process: the PR
  changed **both** together and they must not be reverted independently.

### Analytical validation

The Bofang expression implemented in the C++ process matches the analytical
reproduction exactly (abs error 0) for A, B2 and C; variant B fails because the
assignment is skipped during initialization. See `analytical_verification.py`.

### Classification

**No observable difference** between the pre-#13472 and post-#13472 lifecycle
(faithfully reproduced). The observed difference in variant B is an artifact of
a partial/incorrect revert (dead lifecycle), not of the merged change itself.

### Recommendation

* The current Bofang implementation (`ExecuteBeforeSolutionLoop` +
  `constrained`) is safe for this thermomechanical case: it yields bit-identical
  results to the pre-PR lifecycle and to the historical Kratos.
* Execution should **not** return to `ExecuteInitialize()`. If anything, a
  future refactor should make the analysis robust to the process lifecycle
  method (e.g. drive `ExecuteInitialize()` as well), but the merged behaviour is
  correct.
* A focused regression test is retained (see
  `applications/DamApplication/tests/test_bofang_initialization.py`) that locks
  the physically meaningful invariant: after initialization, the reservoir
  nodes carry the analytical Bofang temperature. This guards the Bofang formula
  and against a future dead-lifecycle regression.

