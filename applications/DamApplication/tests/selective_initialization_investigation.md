# Selective initialization of a new Dam analysis

## Outcome and scope

A fresh `DamAnalysis` can reproduce an equilibrated, nonuniform Q4 gravity
stress field with **zero initial displacement**, using standard constitutive
`InitialState`. It also reproduces the subsequent horizontal-load and thermal
increments of an uninterrupted calculation. This is selective initialization,
not restart continuation.

The working prototype uses two independent Models, core `RestartUtility`, a
small **test-only** C++ bridge, and ordinary process callbacks. No destination
restart path, solver initialization changes, new serializer, nodal stress
mapping, or displacement transfer are needed. Existing reference-temperature
functionality works without a new process.

The work is on `dam/initial-state-transfer` in
`/tmp/kratos-selective-initial-state`. The original checkout remains on
`dam/enable-standard-restart`, with all its uncommitted changes preserved.
Nothing has been committed. The earlier investigation and benchmarks remain.

## Architecture and minimum data

1. Solve Stage A using normal DamAnalysis and standard body acceleration. Save
   a checkpoint through core `save_restart_process`.
2. Create the destination through its normal **mdpa** import and initialize its
   own solver, materials, original geometry, buffer, time, constraints and loads.
3. In an ordinary process's `ExecuteBeforeSolutionLoop`, load the checkpoint
   into a separate source ModelPart using core `RestartUtility`. Do not create
   or initialize a source solver. Keep that Model alive during extraction.
4. Validate same node/element IDs, original coordinates, ordered connectivity,
   SMA element type, integration method/order, reviewed law types and material
   constants. Extract `CAUCHY_STRESS_VECTOR` separately at every source point.
5. Obtain the destination's initialized law pointers using C++
   `CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, ...)`; assign a distinct
   `InitialState(zero_strain, point_stress, identity_F)` to every point.
6. Release the source Model after the copy. Do not repeat assignment during
   `ExecuteInitializeSolutionStep`. Let normal initial GiD output and the new
   analysis proceed.

The **transferred payload** is only three stress components per integration
point, associated with element ID and integration-point order. Matching needs
node IDs/coordinates, connectivity, element/integration identity, law identity
and elastic constants (including thickness, density and thermal expansion where
applicable). These are validation metadata, not destination assignments.

The source restart itself still contains a complete serialized ModelPart:
`RestartUtility` does not selectively deserialize fields. Recovering stresses
requires the source displacement, source temperature/reference temperature and
source law state to be available in that isolated model. None of these are
copied to the destination. No source solve, initialization, load assignment or
state reset is performed. Snapshot comparisons verify unchanged source stresses,
strains, nodal fields, force vectors, time and fixities before/after extraction.
Core loading necessarily sets the **source** `IS_RESTARTED`/`LOAD_RESTART` flags;
this is not modification of the checkpoint or of its mechanical state.

Only requested reference-temperature data would cross the boundary in a thermal
reference-field use case. If already available as node-ID/value data, the
existing reference-temperature process is sufficient; restart is unnecessary.
Current `TEMPERATURE` always belongs to the destination thermal problem.

## Why a small C++ bridge is needed

The existing benchmark showed that assigning the first stress to every Q4
point with `SetInitialStateProcess` causes stress error 4.018, free residual
2.319 and artificial displacement 3.257e-4. A single state per element cannot
represent this field. The C++ interfaces already support distinct point states,
but Python does not expose the law-pointer vector or `SetInitialState`.

`selective_initial_state_probe/bridge.cpp` therefore demonstrates only the
missing connection. It is a separately compiled investigation module, **not a
registered production utility or supported input API**. It uses existing C++
element and law interfaces and owns no serialization, mesh search, interpolation
or restart logic. It validates all elements before assigning any states.

The prototype explicitly rejects incompatible meshes/materials/integration
rules/laws, moved meshes, nonzero destination displacement, an already assigned
initial state, and a thermoelastic destination initially having `T != T_ref`.
The latter restriction makes the initial thermal datum unambiguous. Validated
scope is two-dimensional, fixed-mesh, four-point Q4 SMA small-displacement
plane strain, using `LinearElasticPlaneStrain2DLaw` or Dam
`ThermalLinearElastic2DPlaneStrain`. Other configurations are not enabled.

One independent production prerequisite is included: registration of
`USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN`. Dam writes this variable to
ProcessInfo; an unregistered variable made a real Dam checkpoint unsafe to
load. The existing rebuilt Dam library with this registration was used for the
runs. Its relevant element/law sources match this branch. No conventional
restart lifecycle implementation was merged into this branch.

## Numerical validation through normal DamAnalysis

`test_dam_selective_initialization.py` builds the existing 1 by 2 body using two
unit-square Q4 elements, eight total points, E=10000, nu=0.2, density=2 and
unit thickness. Both base nodes have zero X/Y displacement. Gravity is body
acceleration (0,-10,0), giving total vertical load -40 and reaction +40.
Stage A's maximum displacement is 0.003770008285.

The reference runs gravity at t=1 and t=2, then horizontal acceleration +2 at
t=3. Its t=1 checkpoint is loaded only into the source Model. The fresh
mdpa destination starts at t=0 with zero displacement, reaction, velocity,
acceleration and previous-buffer displacement. Its STEP is 0 and IS_RESTARTED
is false. The test deliberately stamps source checkpoint STEP=41, since this
branch's Dam loop does not increment STEP; the stamp verifies that source
metadata cannot leak into the destination. This instrumentation does not change
mechanics or implement restart step tracking.

The destination retains gravity at t=0, t=1 and t=2 and applies its own
horizontal acceleration at t=3. All comparisons of displacement use
`u_reference(t) - u_reference(gravity)`. Reactions are computed by the new solve;
they start at zero and are not prescribed or transferred.

A second experiment uses the real thermo-mechanical Dam solver with two
`EulerianConvDiff2D4N` thermal elements. The source initially has
`T = T_ref = 20 + 2*y`; the destination independently prescribes
`T = T_ref = 70 + 2*y`. At t=3 both analyses increase only their top/bottom
Dirichlet temperatures by 5. The middle row is free and the thermal solver
computes it: 24 in the reference, 74 in the destination. Nodal conductivity=1,
density=2 and specific heat=1 are supplied, as required by the thermal elements.
Material properties alone do not populate these historical nodal quantities.

Maximum absolute errors over the new initial state and subsequent solves:

| Measure | Gravity + horizontal increment | Gravity + temperature increment |
|---|---:|---:|
| Initial Gauss-point stress difference | 0 | 0 |
| Initial internal/external force difference | 0 / 0 | 0 / 0 |
| Initial residual-vector difference | 0 | 0 |
| First-solve displacement magnitude | 6.17e-19 | 1.26e-18 |
| Displacement increment difference over all steps | 8.78e-19 | 1.66e-18 |
| Stress difference over all steps | 1.28e-14 | 8.88e-15 |
| Reaction difference over all steps | 7.22e-15 | 7.11e-15 |
| Maximum free-DOF residual | 3.55e-15 | 3.55e-15 |

The final new displacement magnitude is 8.203636e-3 for the horizontal increment
and 8.909102e-5 for the thermal increment. Both responses are nontrivial.
The internal-force record is computed as `f_external - residual`; external
forces use the exact constant-body-force integral on these unit-square Q4s.
The residual itself comes directly from assembled element local RHS vectors.
Support entries are balanced by reactions, not included in the free residual.
Removing gravity in a negative test leaves a free residual greater than 1:
matching stresses alone does not ensure equilibrium with arbitrary new loads.

Numeric comparisons require `abs(error) <= 1e-10 + 1e-10*abs(reference)`.
Free residuals must be below 1e-10. Initial zero fields are checked exactly.
The actual initial GiD frame is parsed: all eight stresses match Stage A within
5e-5 absolute (ASCII rounding), and all six nodal displacement vectors are zero.
Later steps prove the states survive the first solve without being reassigned.

## Existing reference-temperature functionality and constitutive convention

The separate `test_existing_reference_temperature_process` needs no bridge or
restart. It uses the existing `impose_nodal_reference_temperature_process` with
a node-ID/value file assigning `T_ref = 20 + 2*y`, and destination current
`T = 30 + 2*y`. All mechanical DOFs are restrained for an analytical check;
thermal top/bottom nodes are prescribed and the middle row is solved. Neither
previous displacement nor previous reaction data exists.

The Dam plane-strain law interpolates both nodal fields at each point and uses
`epsilon_thermal = alpha*(1+nu)*(T-T_ref)*[1,1,0]`.
For alpha=1e-5 and delta T=10 this is `[0.00012,0.00012,0]`; the restrained stress
is `[-5/3,-5/3,0]`. All eight points match these values at t=0 and all three
steps to the stated tolerance. Reference values persist and reactions are
computed on the first solve. This test is registered in the normal small/all
suites; no compiled investigation module is needed for it.

The existing reference process assigns values in `ExecuteBeforeSolutionLoop`
and every `ExecuteInitializeSolutionStep`. Its `interval` setting does not gate
assignment: the regression deliberately uses [0,0] and confirms static reference
values persist. The file's first column is node ID, not time or coordinate.
It changes the reference field, not current temperature or its fixity.

The reviewed elastic implementations subtract `InitialState` strain and add its
stress. With the Dam thermal subtraction, the response is

`stress = C * (strain - initial_strain - thermal_strain) + initial_stress`.

Our selective gravity cases set initial_strain=0 and T=T_ref at initialization.
Thus at zero displacement the exact source stress contributes to internal forces.
Subsequent heating contributes only the new thermal strain. The 50-degree shift
of both destination temperature fields also demonstrates that absolute source
temperatures are not imported and only the intended difference matters.

If preserving the exact imported total stress at t=0 while `T != T_ref` is later
supported, the necessary offset at zero displacement is
`initial_strain = -thermal_strain(t0)`. Alternatively, leaving the offset zero
intentionally superposes an immediate thermal stress and may require a new
equilibrium solve. These are different semantics; the prototype rejects the
ambiguous case rather than silently choosing. Transferring displacements would
require a different offset and is outside this stress-only interface.

## Smallest production proposal

After review of this evidence, promote a narrow same-mesh stress assignment
operation and thin Python binding, with the checks demonstrated by the probe.
A short Python Process should own source `RestartUtility` loading, target
selection, once-only invocation in `ExecuteBeforeSolutionLoop`, and source
lifetime. Keep stress assignment in C++; do not expose mutable source laws to
Python or clone them into the destination. Validate the destination analysis
configuration as quasistatic, fixed mesh, local elasticity before enabling it.

Keep source and destination mechanical model-part names explicit. Allow only
stress transfer in the first version, require zero initial displacement and,
for thermoelastic import, initially equal T and T_ref. Use the existing
reference-temperature process separately, ordered before stress assignment.
Validate the free-DOF initial equilibrium against destination loads and
constraints; incompatible support/loading cannot be cured by a data copy.
The first implementation should fail clearly when preservation of equilibrium
is requested but the normalized residual exceeds its specified tolerance.

No changes are needed to DamAnalysis ordering, destination model import, buffer
filling, time advancement, output scheduling, thermal solver initialization or
restart continuation. No SelfWeight, UP/P, interface, damage or nonlocal changes
are proposed. Reuse standard body-force processes, core restart/save processes,
core `InitialState`, SMA elements, existing Dam thermoelastic laws and the
existing reference-temperature process. Avoid new nodal stress storage,
Gauss-to-node extrapolation, mapping searches, general transfer infrastructure
and custom checkpoint formats. Existing previous-results mapping remains intact.

## Reproduction and evidence

The focused bridge is intentionally outside normal application CMake targets.
From this worktree, using this machine's existing Python 3.14 Kratos build:

```bash
g++ -shared -fPIC -O0 -std=c++20 -fopenmp \
  -DBOOST_UBLAS_MOVE_SEMANTICS -DKRATOS_PYTHON -DKRATOS_SMP_OPENMP \
  -I/usr/include/python3.14 -Ikratos -Iexternal_libraries -I/usr/include/eigen3 \
  -Iapplications/DamApplication -Iapplications/PoromechanicsApplication \
  -Iapplications/StructuralMechanicsApplication -Iapplications/ConstitutiveLawsApplication \
  applications/DamApplication/tests/selective_initial_state_probe/bridge.cpp \
  -L/home/cimne/Kratos/bin/Release/libs \
  -lKratosCore -lKratosStructuralMechanicsCore -lKratosConstitutiveLawsCore \
  -L/home/cimne/Kratos/build/Release/applications/DamApplication -lKratosDamCore \
  -o /tmp/dam_selective_state_probe.so

PYTHONPATH=/tmp:/home/cimne/Kratos/bin/Release \
LD_LIBRARY_PATH=/home/cimne/Kratos/build/Release/applications/DamApplication:/home/cimne/Kratos/bin/Release/libs \
OMP_NUM_THREADS=1 DAM_SELECTIVE_RESULTS=/tmp/dam_selective_results \
python3 applications/DamApplication/tests/test_dam_selective_initialization.py
```

The child analyses explicitly use this worktree's Dam scripts. Rebuild the Dam
core with the variable registration first if using a different build; the probe
must use ABI-compatible libraries. The test fails clearly if the bridge is not
built, rather than silently skipping the selective validation.

Validation on 2026-10-30: all four focused tests passed, including the rejection
subcases; the Dam small suite passed **141 tests** using an isolated Python
runtime with this worktree's Dam scripts. Full source/destination arrays are
recorded in `/tmp/dam_selective_results/horizontal.json` and `thermal.json`;
logs are `/tmp/dam_selective_test.log` and `/tmp/dam_selective_small.log`.
Each run uses temporary source/output directories and leaves input files alone.

The original regression input was read from exactly
`/home/cimne/Simulaciones/00_DEBBUGING_TESTS/DamApplication/TestTranferUtility2D.gid/`
and executed in `/tmp/dam_selective_original_regression`. Its `.post.res` is
byte-identical to `/tmp/dam_standard_restart_regression_final` from the preceding
work. No original case files were changed. This regression is separate from the
controlled selective-initialization benchmark.
