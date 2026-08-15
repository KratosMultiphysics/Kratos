# 2025 DamApplication process-migration audit (PR #13472)

Audit of every DamApplication file changed by PR #13472 ('Dam/assign scalar variable process', merge `359d1ef414`, 2025-06-02, pre-PR parent `ece5cfed146`).

Companion CSV: `2025_process_audit.csv` (same columns).

## Method

- Changed-file inventory reconstructed from git (see below).

- Every C++ process: lifecycle callbacks compared historical vs current (`git show ece5cfe:<file>` vs current); bodies verified identical except the method-name rename (plus `is_fixed`->`constrained` for 7 processes).

- Every Python wrapper: `ExecuteInitialize`->`ExecuteBeforeSolutionLoop` forwarding compared historical vs current.

- `dam_analysis.py`: both process-call sites changed `ExecuteInitialize`->`ExecuteBeforeSolutionLoop` (selfweight block + main block).

- Consumer trace: for each written variable, located the element/condition/constitutive law that reads it and verified it reads only during `CalculateLocalSystem` / `CalculateMaterialResponse` (solve), never during element `Initialize`, `InitializeMaterial`, scheme/strategy `Initialize` or model import.

- Later-fix reconstruction: `git log 359d1ef..HEAD` on the affected files (found #14617).


## Changed-file inventory (PR #13472)

20 C++ processes, 17 Python wrappers, `dam_analysis.py`, 16 test JSON files.

**C++ processes (20):** `apply_component_table_process.hpp`, `dam_added_mass_condition_process.hpp`, `dam_azenha_heat_source_process.hpp`, `dam_bofang_condition_temperature_process.hpp`, `dam_chemo_mechanical_aging_young_process.hpp`, `dam_fix_temperature_condition_process.hpp`, `dam_grouting_reference_temperature_process.hpp`, `dam_hydro_condition_load_process.hpp`, `dam_input_table_nodal_young_modulus_process.hpp`, `dam_nodal_reference_temperature_process.hpp`, `dam_nodal_young_modulus_process.hpp`, `dam_noorzai_heat_source_process.hpp`, `dam_random_fields_variable_process.hpp`, `dam_reservoir_constant_temperature_process.hpp`, `dam_reservoir_monitoring_temperature_process.hpp`, `dam_t_sol_air_heat_flux_process.hpp`, `dam_temperature_by_device_process.hpp`, `dam_uplift_circular_condition_load_process.hpp`, `dam_uplift_condition_load_process.hpp`, `dam_westergaard_condition_load_process.hpp`

**Python wrappers (17):** `apply_constraint_vector_dam_table_process.py`, `apply_load_vector_dam_process.py`, `apply_load_vector_dam_table_process.py`, `dam_analysis.py`, `impose_2d_random_fields_variable_process.py`, `impose_3d_random_fields_variable_process.py`, `impose_chemo_mechanical_aging_process.py`, `impose_face_heat_flux_process.py`, `impose_grouting_reference_temperature_process.py`, `impose_heat_source_process.py`, `impose_input_table_nodal_young_modulus_process.py`, `impose_nodal_reference_temperature_process.py`, `impose_nodal_young_modulus_process.py`, `impose_reservoir_temperature_condition_process.py`, `impose_thermal_parameters_scalar_value_process.py`, `impose_uniform_temperature_process.py`, `impose_water_loads_condition_process.py`

**Tests (16 JSON):** all under `applications/DamApplication/tests/` (`construction`, `joint_*` parameters/results).

## Change-type summary (C++ processes)

| Change type | Processes |
|---|---|

| `is_fixed`->`constrained` AND `ExecuteInitialize`->`ExecuteBeforeSolutionLoop` | apply_component_table, dam_bofang, dam_fix_temperature, dam_nodal_young_modulus, dam_reservoir_constant, dam_reservoir_monitoring |

| `ExecuteInitialize`->`ExecuteBeforeSolutionLoop` only | dam_added_mass, dam_azenha, dam_chemo_mechanical_aging_young, dam_grouting_reference_temperature, dam_hydro, dam_input_table_nodal_young_modulus, dam_nodal_reference_temperature, dam_random_fields, dam_t_sol_air, dam_uplift_circular, dam_uplift, dam_westergaard, dam_noorzai (also renames an internal `this->ExecuteInitialize()` call) |

| `is_fixed`->`constrained` only (callback NOT renamed) | dam_temperature_by_device (init callback unchanged: only `ExecuteInitializeSolutionStep`) |


## Lifecycle call-chain analysis

`dam_analysis.py` drives every process via the process factory. Historical vs current:

```
historical                 current
solver.AddVariables/Import   (same)
processes constructed        (same)
for p: p.ExecuteInitialize()  for p: p.ExecuteBeforeSolutionLoop()   <- pre solver.Initialize
buffer fill                  (same)
gid_output.Initialize        (same)
solver.Initialize()          solver.Initialize()
for p: p.ExecuteBeforeSolutionLoop()  (historical wrappers: no-op,
                                       only ExecuteInitialize existed)
                                       current: re-applies identical values)
time loop: p.ExecuteInitializeSolutionStep()  (both, per step)
```

Every C++ process (except `dam_temperature_by_device`) overrides the init callback, and every wrapper forwards it (historical `ExecuteInitialize`, current `ExecuteBeforeSolutionLoop`). Therefore, for all 19 renamed processes, the value is assigned **before `solver.Initialize()`** in both versions (classification **L0 - lifecycle equivalent**); the current version additionally re-applies the identical value after `solver.Initialize()`.

Because the consumers read the assigned variables only during the solve (verified per family below), the lifecycle rename cannot change the numerical result (same argument proven experimentally for Bofang: A vs B2 vs C bit-identical).


## Assignment-semantics analysis (`ApplyConstantScalarValueProcess` -> `AssignScalarVariableProcess`)

| Aspect | ApplyConstantScalarValueProcess (historical) | AssignScalarVariableProcess (current) |
|---|---|---|

| default fixity | `is_fixed: false` | `constrained: true` |

| Fix behaviour | `Node::Fix` (creates a DOF if missing, then fixes) | `VariableUtils.ApplyFixity` -> `pGetDof->FixDof` (**throws** if the variable has no DOF) |

| Free behaviour | none (applied once) | frees in `ExecuteFinalizeSolutionStep` if it was fixed |

| interval | none (always applied) | `IntervalUtility`; default [0,1e30]; Dam convention [0,0] = active at t=0 only |

| per-step | no re-apply (applied once at `ExecuteInitialize`) | re-applies at every `ExecuteInitializeSolutionStep` inside interval |

| model/modelpart ctor | both | Model only |

| tables | Dam used `ApplyComponentTableProcessDam` for tables | same for Dam; AssignScalarVariableProcess also supports csv tables |

**Consequences**

- For Dam processes that pass `constrained` explicitly, behaviour is equivalent.

- The **default difference** (`is_fixed=false` vs `constrained=true`) is the root cause of the vector-load regression (see below).

- For non-DOF variables with `constrained=true`, historical code silently created a DOF and fixed it (`Node::Fix`); current code **throws** (`ApplyFixity`). Both are wrong for non-DOF quantities, but only the current one surfaces as a hard error.


## Known later regression (vector loads)

- **PR #13472** switched `apply_load_vector_dam_process.py` and `apply_load_vector_dam_table_process.py` from `ApplyConstantScalarValueProcess` (default `is_fixed=false`) to `AssignScalarVariableProcess` without setting `constrained` (default `true`).

- Effect: `ApplyFixity(POINT_LOAD_X,...)` runs on a variable that is not a DOF in the Dam solvers (`AddDofs` only creates DISPLACEMENT/TEMPERATURE). Empirically confirmed on current master: it throws

  `Trying to fix/free dof of variable POINT_LOAD_X but this dof does not exist`.

- **Fix #14617** (2026-07-30, `95e2345ad1a`) added explicit `constrained=false` to all components of both load processes and added `test_apply_load_vector_dam_processes.py` (values + no DOFs).

- Only #14617 touched those files since #13472: the breakage was live for ~13 months and is not covered by any older test.

- Remaining analogous uses audited: every current `AssignScalarVariableProcess` use in DamApplication either passes `constrained` explicitly or forwards it from the settings. No remaining instance of the unguarded default. The `apply_constraint_vector_dam_table_process` correctly forwards `constrained` (a constraint process).


## Newly discovered issue (impose_nodal_young_modulus_process)

- Reproduced on current master: importing `KratosMultiphysics.DamApplication.impose_nodal_young_modulus_process` raises

  `NameError: name 'Process' is not defined`.

- Root cause: PR #13472 removed `from KratosMultiphysics import *` from that wrapper (and updated `Process.__init__` to `KratosMultiphysics.Process.__init__`) but left the class base as bare `Process`, which is now undefined.

- Impact: the wrapper (and therefore the nodal-Young-modulus feature) is completely unusable since #13472. No existing test exercises it (confirmed by searching all DamApplication tests/ProjectParameters).

- Regression demonstration: `tests/test_dam_process_lifecycle.py::test_nodal_young_modulus_wrapper_is_importable` (canary pinning the broken state).

- Per the audit decision rules the fix is **not applied** in this branch; it is reported separately (one line: `class ImposeNodalYoungModulusProcess(KratosMultiphysics.Process):`).


## Full audit table

_Risk: R0 negligible / R1 low / R2 medium / R3 high (see decision rules)._

| Process | Python wrapper | Physical role | Variable(s) written | Entity | Historical callback | Current callback | Lifecycle class (L0/L1/L2/L3) | Historical assignment mechanism | Current assignment mechanism | is_fixed/constrained semantics | Called before solver.Initialize? legacy | Called before solver.Initialize? current | Updated each solution step? | Consumer of value | Initialization consumer? | Existing test coverage | Later fixes | Risk | Recommended action |

|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

| DamBofangConditionTemperatureProcess | impose_reservoir_temperature_condition_process.py (BOFANGTEMPERATURE) | thermal boundary condition (reservoir temperature profile) | TEMPERATURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | default false; wrapper FORCES true (AddEmptyValue+SetBool) -> Fix(TEMPERATURE) on water nodes + Free each step (same pre/post) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | thermal elements (EulerianConvDiff/ConvDiff2D) in CalculateLocalSystem and thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | test_bofang_initialization.py (this branch): analytical + lifecycle + production-wrapper chain; full thermomechanical experiment (A/B2/C) | none | R0 | keep tests; proven numerically identical pre/post #13472 (A vs B2 vs C) |

| DamReservoirConstantTemperatureProcess | impose_reservoir_temperature_condition_process.py (CONSTANTRESERVOIRTEMPERATURE) | thermal boundary condition (constant reservoir temperature) | TEMPERATURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | default false; wrapper FORCES true -> Fix(TEMPERATURE) + Free each step (same pre/post) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | thermal elements (EulerianConvDiff/ConvDiff2D) in CalculateLocalSystem and thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | construction test (impose_reservoir_temperature_condition_process) | none | R0 | none |

| DamReservoirMonitoringTemperatureProcess | impose_reservoir_temperature_condition_process.py (MONITORINGRESERVOIRTEMPERATURE) | thermal boundary condition (measured reservoir temperature) | TEMPERATURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | default false; wrapper FORCES true -> Fix(TEMPERATURE) + Free each step (same pre/post) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | thermal elements (EulerianConvDiff/ConvDiff2D) in CalculateLocalSystem and thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | none | none | R1 | optional small process-level test (assigned value + fixity); lifecycle-equivalent by source |

| DamFixTemperatureConditionProcess | impose_uniform_temperature_process.py; impose_face_heat_flux_process.py (T uniform) | thermal boundary condition (fixed/uniform temperature) | TEMPERATURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | default false; wrapper passes settings value; Fix + Free each step (same pre/post) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | thermal elements (EulerianConvDiff/ConvDiff2D) in CalculateLocalSystem and thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | construction test (impose_uniform_temperature_process) | none | R0 | none |

| DamAzenhaHeatSourceProcess (DamAzenhaHeatFluxProcess) | impose_heat_source_process.py (AzenhaHeatFlux...) | thermal source (hydration heat, Azenha) | HEAT_FLUX, ALPHA_HEAT_SOURCE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep; tracks ALPHA_HEAT_SOURCE history) | thermal elements (ConvDiff2D) volume source in CalculateLocalSystem (solve) | No (only during solve) | none | none | R1 | optional focused unit test; lifecycle-equivalent by source (bodies identical) |

| DamNoorzaiHeatSourceProcess (DamNoorzaiHeatFluxProcess) | impose_heat_source_process.py (NoorzaiHeatFlux...) | thermal source (Noorzai) | HEAT_FLUX | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep; internal call renamed consistently) | thermal elements (ConvDiff2D) volume source in CalculateLocalSystem (solve) | No (only during solve) | none | none | R1 | optional focused unit test; lifecycle-equivalent by source |

| DamTSolAirHeatFluxProcess (DamTSolAirHeatFluxProcess) | impose_face_heat_flux_process.py | thermal boundary condition (T-sol-air ambient heat flux) | FACE_HEAT_FLUX | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | thermal face condition (ConvectionDiffusion) reads nodal FACE_HEAT_FLUX during CalculateLocalSystem (solve); not during Initialize | No (only during solve) | construction test (impose_face_heat_flux_process) | none | R0 | none |

| DamNodalReferenceTemperatureProcess | impose_nodal_reference_temperature_process.py | thermomechanical reference quantity (nodal reference temperature) | NODAL_REFERENCE_TEMPERATURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test (value present before first solve) |

| DamGroutingReferenceTemperatureProcess | impose_grouting_reference_temperature_process.py | thermomechanical reference quantity (grouting) | NODAL_REFERENCE_TEMPERATURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | NO (set at init to mInitialValue; updated only at grouting time in ExecuteFinalizeSolutionStep) | thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test; lifecycle-equivalent by source |

| DamTemperatureByDeviceProcess | impose_temperature_by_device_process.py (temperature_by_device_list) | thermal boundary condition (device-recorded temperature) | TEMPERATURE | nodes | ExecuteInitializeSolutionStep (no init callback; unchanged) | ExecuteInitializeSolutionStep (no init callback; unchanged) | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | N/A (per-step assignment only) | N/A (per-step assignment only) | YES (only at ExecuteInitializeSolutionStep) | thermal elements (EulerianConvDiff/ConvDiff2D) in CalculateLocalSystem and thermomechanical constitutive law (thermal strain) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test |

| DamNodalYoungModulusProcess | impose_nodal_young_modulus_process.py | material parameter (nodal Young modulus) | NODAL_YOUNG_MODULUS | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | default false; wrapper FORCES true -> Fix(NODAL_YOUNG_MODULUS) (Node::Fix creates a DOF; harmless, identical pre/post) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | *Nodal constitutive laws (LinearElastic3DLawNodal / ThermalLinearElastic3DLawNodal) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test; note Fix() on a non-DOF creates a DOF via Node::Fix (no equation impact; identical pre/post) |

| DamInputTableNodalYoungModulusProcess | impose_input_table_nodal_young_modulus_process.py | material parameter (table nodal Young modulus) | NODAL_YOUNG_MODULUS | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep, table by node id) | *Nodal constitutive laws (LinearElastic3DLawNodal / ThermalLinearElastic3DLawNodal) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test |

| DamChemoMechanicalAgingYoungProcess | impose_chemo_mechanical_aging_process.py | material parameter (aging Young modulus) | NODAL_YOUNG_MODULUS | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep; recomputed from time) | *Nodal constitutive laws (LinearElastic3DLawNodal / ThermalLinearElastic3DLawNodal) in CalculateMaterialResponse (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test; note log(time) -> -inf at t=0 (NaN transient overwritten each step; identical pre/post) |

| DamRandomFieldsVariableProcess | impose_2d_random_fields_variable_process.py / impose_3d_random_fields_variable_process.py | random field (material/thermal property) | user-specified (e.g. NODAL_YOUNG_MODULUS, CONDUCTIVITY) | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep, table by node id) | depends on the variable (material or thermal consumers, solve-time) | No (only during solve) | none | none | R1 | optional focused unit test |

| DamHydroConditionLoadProcess | impose_water_loads_condition_process.py (HydroLinePressure2D/HydroSurfacePressure3D) | static load (hydrostatic pressure) | POSITIVE_FACE_PRESSURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | structural surface/line load conditions read nodal POSITIVE_FACE_PRESSURE during CalculateLocalSystem (solve); not during Initialize | No (only during solve) | construction test (impose_water_loads_condition_process) | none | R0 | none |

| DamUpliftConditionLoadProcess | impose_water_loads_condition_process.py (StraightUplift...) | uplift/hydro load | POSITIVE_FACE_PRESSURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | structural surface/line load conditions read nodal POSITIVE_FACE_PRESSURE during CalculateLocalSystem (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test |

| DamUpliftCircularConditionLoadProcess | impose_water_loads_condition_process.py (CircularUplift...) | uplift/hydro load (circular) | POSITIVE_FACE_PRESSURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | structural surface/line load conditions read nodal POSITIVE_FACE_PRESSURE during CalculateLocalSystem (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test |

| DamWestergaardConditionLoadProcess | impose_water_loads_condition_process.py (HydroDynamic...) | dynamic hydrodynamic pressure (Westergaard) | POSITIVE_FACE_PRESSURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | structural surface/line load conditions read nodal POSITIVE_FACE_PRESSURE during CalculateLocalSystem (solve); not during Initialize | No (only during solve) | none | none | R1 | optional focused unit test |

| DamAddedMassConditionProcess | none (instantiated directly; pybind DamAddedMassConditionProcess) | dynamic/added-mass contribution | ADDED_MASS | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | unchanged (parameter rename only) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep) | Dam AddedMassCondition reads nodal ADDED_MASS during CalculateLocalSystem (mass matrix assembly, solve); not during Initialize | No (only during solve) | none | none | R1 | treat separately for dynamic runs; lifecycle-equivalent by source (ADDED_MASS consumed only in CalculateLocalSystem of the added-mass condition during the solve); add a small process-level test if used |

| ApplyComponentTableProcessDam (apply_component_table_process.hpp) | apply_constraint_vector_dam_table_process.py / apply_load_vector_dam_table_process.py (table branches) | generic table assignment (constraint/load components) | user-specified component (DISPLACEMENT_X, POINT_LOAD_X, ...) | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (effective assignment pre solver.Initialize in both versions; consumers read only during solve) | FastGetSolutionStepValue(var) = value (same body, method renamed only) | FastGetSolutionStepValue(var) = value (identical body) | default false; Fix only if constrained (same pre/post) | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (ExecuteInitializeSolutionStep reads table by time) | depends on variable (DOF or load; solve-time) | No (only during solve) | test_apply_load_vector_dam_processes.py (table path) | none | R1 | keep existing test |

| ImposeNodalYoungModulusProcess (wrapper) | impose_nodal_young_modulus_process.py | material parameter (nodal Young modulus) | NODAL_YOUNG_MODULUS | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (wrapper forwards the renamed callback; value assigned pre solver.Initialize in both versions) | DamNodalYoungModulusProcess via wrapper (is_fixed forced true) | DamNodalYoungModulusProcess via wrapper (constrained forced true) | wrapper forces constrained=true (same pre/post) | YES (pre-solver-init via ExecuteInitialize) | N/A - wrapper is unimportable | would be YES (per-step) once importable | nodal constitutive laws (solve) | No | none | none | R3 - NEWLY DISCOVERED: wrapper cannot even be imported | fix the class base (Process -> KratosMultiphysics.Process); then enable the functional harness test; reported separately (not applied in this branch, per decision rules) |

| ApplyConstraintVectorDamTableProcess (wrapper) | apply_constraint_vector_dam_table_process.py | constraint (prescribed displacement) | DISPLACEMENT_X / _Y / _Z | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (wrapper forwards the renamed callback; value assigned pre solver.Initialize in both versions) | ApplyConstantScalarValueProcess(model_part, params), is_fixed passed from JSON (default false) | AssignScalarVariableProcess(Model, params), constrained passed from JSON; ApplyComponentTableProcessDam for tables | explicitly passed from JSON (intended: true for constraints); identical values pre/post | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (AssignScalarVariableProcess re-applies each step) | DISPLACEMENT is the primary DOF; consumed by the solver (fixed value on constrained DOFs) | No (DOF values read at solve; fixity set before first SetUpDofSet in both versions) | joint tests (all 8), construction test, bofang case | none needed (constrained passed explicitly) | R0 | none (covered by existing tests) |

| ApplyLoadVectorDamProcess (wrapper) | apply_load_vector_dam_process.py | static load (point-load vector) | POINT_LOAD_X / _Y / _Z | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (wrapper forwards the renamed callback; value assigned pre solver.Initialize in both versions) | ApplyConstantScalarValueProcess(model_part, params); is_fixed NOT set -> default false -> no fix, value set once | AssignScalarVariableProcess(Model, params); post-PR constrained NOT set -> default TRUE -> ApplyFixity on POINT_LOAD_* (no DOF) -> KRATOS_ERROR crash; FIXED by #14617 (explicit constrained=false) | default differs! ApplyConstantScalarValueProcess default is_fixed=false; AssignScalarVariableProcess default constrained=true. #14617 sets false explicitly. | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (AssignScalarVariableProcess re-applies each step) | structural point-load condition reads POINT_LOAD during CalculateLocalSystem (solve) | No | test_apply_load_vector_dam_processes.py (added by #14617) | #14617 (2026-07-30): explicit constrained=false for all components | R2 historically (real breakage between #13472 and #14617); current R0/R1 (fixed + tested) | keep regression test; documented in audit |

| ApplyLoadVectorDamTableProcess (wrapper) | apply_load_vector_dam_table_process.py | static load (table point-load vector) | POINT_LOAD_X / _Y / _Z | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (wrapper forwards the renamed callback; value assigned pre solver.Initialize in both versions) | ApplyConstantScalarValueProcess + ApplyComponentTableProcessDam (tables) | AssignScalarVariableProcess (constrained=false after #14617) + ApplyComponentTableProcessDam (tables) | same default-difference as ApplyLoadVectorDamProcess; fixed by #14617 | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES (table read each step) | structural point-load condition, CalculateLocalSystem (solve) | No | test_apply_load_vector_dam_processes.py (added by #14617) | #14617 (explicit constrained=false) | R2 historically; current R0/R1 (fixed + tested) | keep regression test; documented in audit |

| ImposeThermalParametersScalarValueProcess (wrapper) | impose_thermal_parameters_scalar_value_process.py | material/thermal parameter | DENSITY, CONDUCTIVITY, SPECIFIC_HEAT | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (wrapper forwards the renamed callback; value assigned pre solver.Initialize in both versions) | ApplyConstantScalarValueProcess(model_part, params) (is_fixed not set -> false) | AssignScalarVariableProcess(Model, params) with explicit constrained=false | explicit false (not DOFs); equivalent | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES | thermal elements read nodal CONDUCTIVITY/SPECIFIC_HEAT/DENSITY in CalculateLocalSystem (solve) | No | construction test | none | R0 | none |

| ImposeUniformTemperatureProcess (wrapper) | impose_uniform_temperature_process.py | thermal boundary condition (uniform temperature) | TEMPERATURE | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (wrapper forwards the renamed callback; value assigned pre solver.Initialize in both versions) | ApplyConstantScalarValueProcess / DamFixTemperatureConditionProcess | constrained from settings: true -> DamFixTemperatureConditionProcess; false -> AssignScalarVariableProcess(constrained from settings) | passed through; equivalent | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES | thermal elements + thermo CL (solve) | No | construction test | none | R0 | none |

| ImposeFaceHeatFluxProcess (wrapper) | impose_face_heat_flux_process.py | thermal boundary condition (uniform T + T-sol-air flux) | TEMPERATURE (uniform), FACE_HEAT_FLUX (T-sol-air) | nodes | ExecuteInitialize | ExecuteBeforeSolutionLoop | L0 (wrapper forwards the renamed callback; value assigned pre solver.Initialize in both versions) | ApplyConstantScalarValueProcess (is_fixed=false) + DamFixTemperatureConditionProcess + DamTSolAirHeatFluxProcess | AssignScalarVariableProcess (constrained=false) + DamFixTemperatureConditionProcess + DamTSolAirHeatFluxProcess | explicit false for the uniform-T component; equivalent | YES (pre-solver-init via ExecuteInitialize) | YES (pre-solver-init, value-identical; current also re-applies post-solver-init) | YES | thermal elements + thermal face condition (solve) | No | construction test | none | R0 | none |


## Risk summary

| Risk | count |
|---|---|

| R0 | 9 |

| R1 | 15 |

| R2 | 2 |

| R3 | 1 |


- **R0 (9):** Bofang, reservoir_constant, fix_temperature, TSolAir, hydro, apply_constraint_vector (wrapper), impose_thermal_parameters (wrapper), impose_uniform_temperature (wrapper), impose_face_heat_flux (wrapper).

- **R1 (15):** reservoir_monitoring, azenha, noorzai, nodal_reference_temperature, grouting_reference_temperature, temperature_by_device, nodal_young_modulus, input_table_nodal_young_modulus, chemo_mechanical_aging, random_fields, uplift, uplift_circular, westergaard, added_mass, apply_component_table. Lifecycle-equivalent by source but untested.

- **R2 (2, historical only):** apply_load_vector / apply_load_vector_table — real breakage between #13472 and #14617; now fixed and tested.

- **R3 (1):** impose_nodal_young_modulus_process (wrapper) — newly discovered, unimportable since #13472; fix reported separately (see dedicated section).


## Decisions / recommendations

- No production-code change is required by this audit for the lifecycle/assignment semantics: all are equivalent for the merged behaviour.

- The only real semantic regression found (#13472 load fixity) is already fixed by #14617 and has regression tests.

- **One production fix is recommended (reported separately, not applied here):** `impose_nodal_young_modulus_process.py` class base `Process` -> `KratosMultiphysics.Process` (one line).

- Optional focused tests (R1 rows) can reuse the generic process harness below; they are not required for correctness but protect the per-process lifecycle coupling.
