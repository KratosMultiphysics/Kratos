---
title: Processes
keywords: contact, processes, contact search, ALMContactProcess, SearchBaseProcess, ALMFastInit, NormalGapProcess, NormalCheckProcess, AALM, dynamic factor, master slave, ProcessFactoryUtility
tags: [contact, implementation, processes, search, ALM, MPC, mesh tying]
sidebar: contact_structural_mechanics_application
summary: Reference of the sixteen C++ processes and nine Python processes of the ContactStructuralMechanicsApplication, with their purpose, constructors and verbatim default JSON, Python names, the stage of the solution pipeline in which they run, the Execute* life cycle of the Python contact processes and the call sequence between them.
---

> **Sources.** Thesis §4.4 (contact search, pp. 123–137), App. D.4.3.1 (adapted augmented Lagrangian, Algorithm 7, p. 320); code: `custom_processes/*.{h,cpp}`, `custom_python/add_custom_processes_to_python.cpp`, `custom_python/process_factory_utility.{h,cpp}`, `python_scripts/search_base_process.py`, `python_scripts/alm_contact_process.py`, `python_scripts/penalty_contact_process.py`, `python_scripts/explicit_penalty_contact_process.py`, `python_scripts/mpc_contact_process.py`, `python_scripts/mesh_tying_process.py`, `python_scripts/contact_remesh_mmg_process.py`, `python_scripts/basic_mapping_process.py`, `python_scripts/replace_properties_process.py`, `custom_strategies/custom_convergencecriterias/base_mortar_criteria.h`; tests `tests/cpp_tests/processes/*.cpp`, `tests/test_dynamic_search.py`, `tests/test_check_normals_process.py`, `tests/test_process_factory.py`.

A *process* in Kratos is an object with the standard hooks `ExecuteInitialize`, `ExecuteBeforeSolutionLoop`, `ExecuteInitializeSolutionStep`, `ExecuteFinalizeSolutionStep`, `ExecuteBeforeOutputStep`, `ExecuteAfterOutputStep` and `ExecuteFinalize`, called by the analysis stage at the corresponding moments of the simulation. The ContactStructuralMechanicsApplication uses processes at two levels:

- **Python processes** (`python_scripts/*_process.py`) are the user-facing entry points listed in `ProjectParameters.json` (`contact_process_list`). They build the contact model parts, set the process-info values and own the C++ search objects.
- **C++ processes** (`custom_processes/`) perform the heavy lifting: the search for contact pairs, the gap mapping, the normal check, the ALM initialization and the penalty adaptation. Most of them are exposed to Python and can also be driven directly (the C++ unit tests and `tests/test_dynamic_search.py` do so).

This page is the *implementation* reference of both layers. The search algorithms themselves (broad phase, OBB/SAT narrow phase, pair filtering, activation) are described in [Search pipeline and bounding volumes](../Contact_Search/Search_Pipeline_And_Bounding_Volumes.html) and [Gap computation](../Contact_Search/Gap_Computation.html); the user-facing JSON of the contact processes is documented key by key in the [Contact process settings reference](../Usage/Contact_Process_Settings_Reference.html). The place of the processes inside the time step is summarized in [Architecture](Architecture.html).

## C++ processes at a glance

All classes live in `custom_processes/` and derive from `Kratos::Process`. The "Python name" column is the name registered in [`add_custom_processes_to_python.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_python/add_custom_processes_to_python.cpp) inside the `KratosMultiphysics.ContactStructuralMechanicsApplication` module (abbreviated `CSMA` below).

| Class | Files | Python name | Constructor(s) exposed | When it runs |
|---|---|---|---|---|
| `BaseContactSearchProcess<TDim,TNumNodes,TNumNodesMaster>` | `base_contact_search_process.{h,cpp}` | (abstract, not exposed) | – | Owns the search life cycle used by all derived search processes |
| `SimpleContactSearchProcess<…>` | `simple_contact_search_process.{h,cpp}` | `SimpleContactSearchProcess2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N` | `(ModelPart&)`, `(ModelPart&, Parameters)`, `(ModelPart&, Parameters, Properties::Pointer)` | Every search (`simple_search: true`) |
| `AdvancedContactSearchProcess<…>` | `advanced_contact_search_process.{h,cpp}` | `AdvancedContactSearchProcess2D2N`, … `3D4N3N` | same three | Every search (`simple_search: false`, default) |
| `ContactSearchWrapperProcess` | `contact_search_wrapper_process.{h,cpp}` | `ContactSearchProcess` | same three | Created by `SearchBaseProcess._create_main_search`; forwards every hook |
| `MPCContactSearchProcess<…>` | `mpc_contact_search_process.{h,cpp}` | `MPCContactSearchProcess2D2N`, … `3D4N3N` | same three | Every search of the MPC formulation |
| `MPCContactSearchWrapperProcess` | `mpc_contact_search_wrapper_process.{h,cpp}` | `MPCContactSearchProcess` | same three | Created by `MPCContactProcess._create_main_search` |
| `NormalGapProcess<TDim,TNumNodes,TNumNodesMaster>` | `normal_gap_process.{h,cpp}` | `NormalGapProcess2D2N`, … `3D4N3N` | `(ModelPart& master, ModelPart& slave)`, `(…, const bool SearchOrientation)` | Inside the search (`check_gap: MappingCheck`) |
| `NormalCheckProcess` | `normal_check_process.{h,cpp}` | `NormalCheckProcess` | `(ModelPart&)`, `(ModelPart&, Parameters)` | Once, in `SearchBaseProcess.ExecuteInitialize` (not on restart) |
| `ALMFastInit` | `alm_fast_init_process.{h,cpp}` | `ALMFastInit` | `(ModelPart&)` | Once, in `_initialize_search_conditions` of the ALM, penalty and MPC processes |
| `ALMVariablesCalculationProcess` | `alm_variables_calculation_process.{h,cpp}` | `ALMVariablesCalculationProcess` | `(ModelPart&)`, `(ModelPart&, Variable<double>&)`, `(ModelPart&, Variable<double>&, Parameters)` | Once, in `_initialize_problem_parameters` (unless `manual_ALM`) and in `MeshTyingProcess.ExecuteInitialize` |
| `AALMAdaptPenaltyValueProcess` | `aalm_adapt_penalty_value_process.{h,cpp}` | (not exposed) | `(ModelPart&)` | Every non-linear iteration from `BaseMortarConvergenceCriteria::PreCriteria` when `ADAPT_PENALTY` |
| `ComputeDynamicFactorProcess` | `compute_dynamic_factor_process.{h,cpp}` | `ComputeDynamicFactorProcess` | `(ModelPart&)` | Every iteration from `PreCriteria` in dynamic problems with `compute_dynamic_factor`; every step from `ExplicitPenaltyContactProcess` |
| `MasterSlaveProcess` | `master_slave_process.{h,cpp}` | `MasterSlaveProcess` | `(ModelPart&)` | Once, when the `Contact` sub-model-part already exists |
| `AssignParentElementConditionsProcess` | `assign_parent_element_conditions_process.{h,cpp}` | `AssignParentElementConditionsProcess` | `(ModelPart& conditions, ModelPart& elements)`, `(Model&, Parameters)` | `MeshTyingProcess` with `consider_static_condensation` |
| `FindIntersectedGeometricalObjectsWithOBBContactSearchProcess` | `find_intersected_geometrical_objects_with_obb_for_contact_search_process.{h,cpp}` | `FindIntersectedGeometricalObjectsWithOBBContactSearchProcess` | `(ModelPart&, ModelPart&)`, `(ModelPart&, ModelPart&, const double)`, `(Model&, Parameters)` | Inside the search when `type_search` is `OctreeWithOBB` |
| `ContactSPRErrorProcess<TDim>` | `contact_spr_error_process.{h,cpp}` | `ContactSPRErrorProcess2D`, `ContactSPRErrorProcess3D` | `(ModelPart&)`, `(ModelPart&, Parameters)` | Before remeshing, created by `ContactRemeshMmgProcess._GenerateErrorProcess` |

Only the processes with a `GetDefaultParameters` method accept a `Parameters` object; `ALMFastInit`, `MasterSlaveProcess`, `ComputeDynamicFactorProcess`, `AALMAdaptPenaltyValueProcess` and `NormalGapProcess` are configured exclusively through the model part (flags, `ProcessInfo` values and nodal variables).

## The search process family

### `BaseContactSearchProcess`

[`base_contact_search_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/base_contact_search_process.h) is the class that actually searches for contact pairs and creates the paired mortar conditions in the `ComputingContact` sub-model-part. It is templated on the dimension, the number of slave nodes and the number of master nodes and explicitly instantiated for `<2,2>`, `<3,3>`, `<3,4>`, `<3,3,4>` and `<3,4,3>`. Its constructor takes the *main* model part (which must contain the `Contact` sub-model-part), a `Parameters` object and an optional `Properties::Pointer` used for the paired conditions. The default parameters are:

```json
{
    "allocation_size"                      : 1000,
    "bucket_size"                          : 4,
    "search_factor"                        : 3.5,
    "type_search"                          : "InRadius",
    "check_gap"                            : "MappingCheck",
    "condition_name"                       : "",
    "final_string"                         : "",
    "inverted_search"                      : false,
    "dynamic_search"                       : false,
    "static_check_movement"                : false,
    "predefined_master_slave"              : true,
    "id_name"                              : "",
    "normal_orientation_threshold"         : 1.0e-1,
    "consider_gap_threshold"               : false,
    "predict_correct_lagrange_multiplier"  : false,
    "pure_slip"                            : false,
    "debug_mode"                           : false,
    "octree_search_parameters" : {
        "bounding_box_factor"             : 0.1,
        "debug_obb"                       : false,
        "OBB_intersection_type"           : "SeparatingAxisTheorem",
        "build_from_bounding_box"         : true
    }
}
```

The meaning of every key is given in the [`search_parameters` reference](../Contact_Search/Search_Pipeline_And_Bounding_Volumes.html#search_parameters-reference). The Python layer fills `condition_name`, `final_string`, `id_name`, `predefined_master_slave` and `pure_slip` automatically (see `_create_search_parameters` below); the user never writes them.

**Enumerations.** The header defines four `enum class` types that are used throughout the search code:

| Enum | Values | Meaning |
|---|---|---|
| `SearchTreeType` | `KdtreeInRadius = 0`, `KdtreeInBox = 1`, `KdtreeInRadiusWithOBB = 2`, `KdtreeInBoxWithOBB = 3`, `OctreeWithOBB = 4`, `Kdop = 5` | Broad-phase structure selected by `type_search` (`ConvertSearchTree`). `Kdop` is not implemented and raises an error |
| `CheckResult` | `Fail = 0`, `AlreadyInTheMap = 1`, `OK = 2` | Result of `CheckGeometricalObject` / `CheckCondition` for a candidate pair |
| `CheckGap` | `NoCheck = 0`, `DirectCheck = 1`, `MappingCheck = 2` | Pair-creation policy selected by `check_gap` (`ConvertCheckGap`), see [`check_gap` modes](../Contact_Search/Gap_Computation.html#check_gap-modes) |
| `TypeSolution` | `NormalContactStress = 0`, `ScalarLagrangeMultiplier = 1`, `VectorLagrangeMultiplier = 2`, `FrictionlessPenaltyMethod = 3`, `FrictionalPenaltyMethod = 4`, `OtherFrictionless = 5`, `OtherFrictional = 6` | Formulation detected in the constructor from the nodal variables (`LAGRANGE_MULTIPLIER_CONTACT_PRESSURE`, `SCALAR_LAGRANGE_MULTIPLIER`, `VECTOR_LAGRANGE_MULTIPLIER`, `WEIGHTED_GAP`); it decides which multiplier is zeroed in `ClearMortarConditions` and initialized in `SetActiveNode` |

**Local flags.** Five `KRATOS_DEFINE_LOCAL_FLAG` flags store the boolean configuration on the process object itself: `INVERTED_SEARCH` (from `inverted_search`, toggled at run time by `InvertSearch()`), `CREATE_AUXILIAR_CONDITIONS` (`true` when `condition_name` is not empty, i.e. paired conditions are created), `MULTIPLE_SEARCHS` (`true` when `id_name` is not empty, i.e. the process works on `ContactSub<id_name>` and `ComputingContactSub<id_name>` instead of on the whole `Contact` model part), `PREDEFINE_MASTER_SLAVE` (from `predefined_master_slave`) and `PURE_SLIP` (from `pure_slip`). The inquiry helpers `IsInvertedSearch`, `IsMultipleSearchs`, `IsPureSlip` and their `IsNot*` counterparts wrap them. Two compile-time constants complete the configuration: `GapThreshold = 2.0e-3` (simple activation, multiplied by `NODAL_H`) and `ZeroTolerance = std::numeric_limits<double>::epsilon()`.

**Public API.** `Execute()` calls the three hooks in sequence; the individual hooks and the methods they call are all exposed to Python:

| Method | Called by | Work |
|---|---|---|
| `ExecuteInitialize()` | `SearchBaseProcess.ExecuteInitialize` | `CheckContactModelParts()` → `CreatePointListMortar()` → `InitializeMortarConditions()` |
| `ExecuteInitializeSolutionStep()` | `SearchBaseProcess.ExecuteInitializeSolutionStep` (when `_compute_search()`) | `ClearMortarConditions()` → `UpdateMortarConditions()` |
| `ExecuteFinalizeSolutionStep()` | `SearchBaseProcess.ExecuteFinalizeSolutionStep` | `ClearMortarConditions()` |
| `CheckContactModelParts()` | `ExecuteInitialize` | Clones the conditions of `ContactSub<N>` flagged `MARKER` (conditions shared by several pairs) so that each pair has its own `INDEX_MAP` |
| `CreatePointListMortar()` / `UpdatePointListMortar()` | `ExecuteInitialize` / `UpdateMortarConditions` | Fill / refresh the kd-tree point list with the destination conditions (masters, or slaves for an inverted search); in the dynamic case the centers are moved with `ContactUtilities::ComputeStepJump` |
| `InitializeMortarConditions()` | `ExecuteInitialize` | Gives every condition of `ContactSub<N>` an empty `INDEX_MAP` |
| `ClearMortarConditions()` | both step hooks | Calls `ResetContactOperators()` and zeroes the multiplier of the inactive nodes (`ClearScalarMortarConditions`, `ClearComponentsMortarConditions` or `ClearALMFrictionlessMortarConditions` according to `TypeSolution`) |
| `ResetContactOperators()` | `ClearMortarConditions` (virtual, overridden by the MPC variant) | Removes from `ComputingContactSub<N>` the paired conditions that are `ACTIVE = false` (or all of them after remeshing, root flagged `MODIFIED`) and their ids from the `INDEX_MAP` |
| `UpdateMortarConditions()` | `ExecuteInitializeSolutionStep` | The actual search: broad phase, narrow phase, `CheckPairing`, `ComputeActiveInactiveNodes` |
| `CheckMortarConditions()` | debugging only (commented out in `ExecuteInitializeSolutionStep`) | Prints the `INDEX_MAP` of every condition and the active nodes |
| `InvertSearch()` | user code | Flips `INVERTED_SEARCH` |

The virtual methods `CheckPairing`, `ComputeActiveInactiveNodes`, `SetActiveNode`, `SetInactiveNode`, `CleanModelPart` and `AddPairing` are the extension points used by the derived classes.

### `SimpleContactSearchProcess`

[`simple_contact_search_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/simple_contact_search_process.h) only overrides `SetActiveNode(Node&, CommonEpsilon, ScaleFactor)`: the node is activated with the base criterion $$g_n \lt 2\cdot 10^{-3}\,h$$ and, if it penetrates, its multiplier is initialized with the penalty guess $$\lambda = \varepsilon A_i g_n / k$$ (see [Gap computation](../Contact_Search/Gap_Computation.html#activation-thresholds)). It is selected with `"simple_search": true` and is the cheapest activation strategy, suited to problems where the initial active set is easy to guess (patch tests, flat interfaces).

### `AdvancedContactSearchProcess`

[`advanced_contact_search_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/advanced_contact_search_process.h) is the default (`"simple_search": false`). It overrides `CheckPairing` (computes `DISTANCE_THRESHOLD` as the maximum of the mean `NODAL_H` of the two sides before mapping the gap with `NormalGapProcess`), `ComputeActiveInactiveNodes` (activation against `ACTIVE_CHECK_FACTOR * DISTANCE_THRESHOLD`, the weighted gap and the `static_check_movement` / `consider_gap_threshold` corrections, with its own `GapThreshold = 2.0e-4`) and adds the Lagrange multiplier prediction/correction machinery: `ComputeLinearRegressionGapPressure` fits a pressure–gap line over the active nodes and `SetActiveNodeWithRegression` applies `Predict*MortarLM` / `Correct*MortarLM` (`Scalar`, `Components`, `ALMFrictionless`, `ALMFrictionlessComponents`, `ALMFrictional` variants) when `predict_correct_lagrange_multiplier` is `true`. The algorithm is given in pseudo-code in [Gap computation](../Contact_Search/Gap_Computation.html#activation-thresholds).

### `ContactSearchWrapperProcess`

[`contact_search_wrapper_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/contact_search_wrapper_process.h) is the class the Python layer instantiates (`CSMA.ContactSearchProcess`). It removes the need to know the template instance: its constructor reads `DOMAIN_SIZE` from the `ProcessInfo`, loops over the conditions of the model part to get the number of nodes of the first `MASTER` (`size_1`) and the first `SLAVE` (`size_2`) condition (or of any condition when `predefined_master_slave` is `false`), reads and removes `simple_search`, and creates the matching `SimpleContactSearchProcess` or `AdvancedContactSearchProcess` instance (`<2,2>`, `<3,3>`, `<3,4>`, `<3,3,4>` for master quadrilaterals / slave triangles, `<3,4,3>` for the opposite). `Execute`, `ExecuteInitialize`, `ExecuteInitializeSolutionStep` and `ExecuteFinalizeSolutionStep` are forwarded to the wrapped process. Its defaults are the base ones plus `simple_search` and the two asymmetric OBB coefficients:

```json
{
    "simple_search"                        : false,
    "allocation_size"                      : 1000,
    "bucket_size"                          : 4,
    "search_factor"                        : 3.5,
    "type_search"                          : "InRadius",
    "check_gap"                            : "MappingCheck",
    "condition_name"                       : "",
    "final_string"                         : "",
    "inverted_search"                      : false,
    "dynamic_search"                       : false,
    "static_check_movement"                : false,
    "predefined_master_slave"              : true,
    "id_name"                              : "",
    "normal_orientation_threshold"         : 1.0e-1,
    "consider_gap_threshold"               : false,
    "predict_correct_lagrange_multiplier"  : false,
    "pure_slip"                            : false,
    "debug_mode"                           : false,
    "octree_search_parameters" : {
        "bounding_box_factor"             : 0.1,
        "debug_obb"                       : false,
        "OBB_intersection_type"           : "SeparatingAxisTheorem",
        "build_from_bounding_box"         : true,
        "lower_bounding_box_coefficient"  : 0.0,
        "higher_bounding_box_coefficient" : 1.0
        }
}
```

### `MPCContactSearchProcess` and `MPCContactSearchWrapperProcess`

[`mpc_contact_search_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/mpc_contact_search_process.h) derives from `BaseContactSearchProcess` and serves the constraint-based formulation described in [Frictional laws and MPC constraint](Frictional_Laws_And_MPC_Constraint.html). Its `AddPairing` override first creates the paired `MPCMortarContactCondition` through the base implementation and then creates a `ContactMasterSlaveConstraint` (id = maximum constraint id + 1), sets it `ACTIVE`, initializes it, adds it to the computing model part and stores it in the condition under `CONSTRAINT_POINTER`; the condition additionally receives the `SLIP` flag (frictional problems, main model part `SLIP`) or the `RIGID` flag (main model part `RIGID`). Because constraints must be cleaned together with the conditions, it also overrides `CheckContactModelParts` (clones `MARKER` constraints), `ResetContactOperators` (flags the constraints of the inactive pairs `TO_ERASE` and removes them with `RemoveMasterSlaveConstraintsFromAllLevels`) and `CleanModelPart`. The wrapper [`mpc_contact_search_wrapper_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/mpc_contact_search_wrapper_process.h) (Python name `MPCContactSearchProcess`) has exactly the same defaults as `ContactSearchWrapperProcess` and the same dispatch logic, but always creates an `MPCContactSearchProcess` (the `simple_search` key is accepted, removed and ignored).

## Gap and normal processes

### `NormalGapProcess`

[`normal_gap_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/normal_gap_process.h) computes the consistent nodal gap `NORMAL_GAP` (thesis Algorithm 8, see [Gap computation](../Contact_Search/Gap_Computation.html#consistent-nodal-gap-normal_gap-thesis-algorithm-8)). It is constructed with the master and slave model parts (`MasterSubModelPart<N>` / `SlaveSubModelPart<N>`, created by `BaseContactSearchProcess::SetOriginDestinationModelParts`) and an optional `SearchOrientation` flag (`true` by default; `false` for an inverted search, in which case the `MASTER`/`SLAVE` node flags are swapped with `SwitchFlagNodes` for the duration of the mapping). `Execute()`:

1. stores the coordinates of the origin nodes in `AUXILIAR_COORDINATES` and zeroes them on the destination side;
2. maps `AUXILIAR_COORDINATES` from master to slave with the core `SimpleMortarMapperProcess` (`MapperType`), forwarding `DISTANCE_THRESHOLD` and `ZERO_TOLERANCE_FACTOR` from the `ProcessInfo` and `CONSIDER_TESSELLATION` from the slave properties, with `remove_isolated_conditions: true` and `update_interface: false`;
3. `ComputeNormalGap` sets, on every slave node whose mapped coordinates are not zero, `NORMAL_GAP` $$= (\mathbf{x}_s - \mathbf{x}_{mapped}) \cdot (-\mathbf{n}_s)$$; nodes on the other side get `NORMAL_GAP = 0`.

It is created and executed inside `BaseContactSearchProcess::ComputeMappedGap` (called from `CheckPairing` when `check_gap` is `MappingCheck`) after `NORMAL_GAP` has been initialized to `1.0e12`; it is not meant to be used as a standalone process in `ProjectParameters.json`, although the Python bindings (`NormalGapProcess2D2N`, …) allow it.

### `NormalCheckProcess`

[`normal_check_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/normal_check_process.h) verifies that the normals of the solid element faces (and of the conditions on them) point outwards and inverts the wrong ones. Defaults:

```json
{
    "length_proportion" : 0.1,
    "check_threshold"   : 5.0e-7
}
```

The algorithm (stage 2 of the [search pipeline](../Contact_Search/Search_Pipeline_And_Bounding_Volumes.html#stage-2--normal-check-normalcheckprocess)) computes the normals with `NormalCalculationUtils`, moves the center of each face by `length_proportion` times its length along the normal and tests `IsInside` on the parent element with `check_threshold`; faces whose offset center falls inside the element are flagged `MARKER` together with their nodes, the conditions sharing those nodes inherit the flag, and `MortarUtilities::InvertNormalForFlag` inverts the geometry of the flagged elements and conditions. Slender elements (beams, shells, membranes) are skipped with an informative message; all `MARKER` flags are reset at the end. `SearchBaseProcess.ExecuteInitialize` runs it once with `length_proportion = normal_check_proportion` unless `IS_RESTARTED` is set. Tests: `tests/test_check_normals_process.py` (`test_check_normals`, `test_check_normals_quads`, `test_check_normals_s_shape`).

## ALM initialization and adaptation processes

### `ALMFastInit`

[`alm_fast_init_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/alm_fast_init_process.h) is constructed with the `Contact` model part and executed once by `_initialize_search_conditions` of `ALMContactProcess`, `PenaltyContactProcess`, `ExplicitPenaltyContactProcess` and `MPCContactProcess`. `Execute()`:

1. renumbers the conditions of the **root** model part consecutively from 1 (several utilities assume ordered ids);
2. reads the frictional flag from the model part (`Is(SLIP)`) and the global penalty from `ProcessInfo[INITIAL_PENALTY]`;
3. on every node that is `SLAVE` (or whose `SLAVE` flag is undefined) sets `WEIGHTED_GAP = 0`, `WEIGHTED_SLIP = 0` (frictional), nodal `INITIAL_PENALTY` $$= \varepsilon$$, `DYNAMIC_FACTOR = 1`, `AUGMENTED_NORMAL_CONTACT_PRESSURE = 0` and `AUGMENTED_TANGENT_CONTACT_PRESSURE = 0` (frictional);
4. sets the condition `NORMAL` to zero (it is recomputed by the search);
5. in frictional problems computes the nodal `FRICTION_COEFFICIENT` as the average over the conditions attached to the node of the `FRICTION_COEFFICIENT` of their properties (a warning is printed for properties without it). `NODAL_AREA` is used as a temporary counter here and is overwritten by the first explicit contribution.

### `ALMVariablesCalculationProcess`

[`alm_variables_calculation_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/alm_variables_calculation_process.h) computes the global penalty $$\varepsilon$$ (`INITIAL_PENALTY`) and scale factor $$k$$ (`SCALE_FACTOR`) of the augmented Lagrangian from the material and the mesh. Constructor: `(ModelPart&, const Variable<double>& rNodalLengthVariable = NODAL_H, Parameters = {})`; the nodal length variable must be a historical variable of the model part. Defaults:

```json
{
    "stiffness_factor"     : 10.0,
    "penalty_scale_factor" : 1.0,
    "compute_scale_factor" : true,
    "compute_penalty"      : true
}
```

`Execute()` loops over the conditions of the contact model part and accumulates, separately for `SLAVE` and `MASTER` conditions (conditions without either flag count on both sides), the domain size, the volume-weighted `YOUNG_MODULUS` of the condition properties and the area-weighted nodal length. With the mean values $$\bar{E}$$ and $$\bar{h}$$ of each side it sets

<p align="center">$$ \varepsilon = \min\left( s \frac{\bar{E}_{slave}}{\bar{h}_{slave}},\; s \frac{\bar{E}_{master}}{\bar{h}_{master}} \right), \qquad k = \min\left( p\, s \frac{\bar{E}_{slave}}{\bar{h}_{slave}},\; p\, s \frac{\bar{E}_{master}}{\bar{h}_{master}} \right) $$</p>

with $$s$$ = `stiffness_factor` and $$p$$ = `penalty_scale_factor` (the *less stiff* side wins; if the slave side has no Young modulus the master value is used for $$\varepsilon$$). `ALMContactProcess._initialize_problem_parameters` calls it with `stiffness_factor` and `penalty_scale_factor` taken from `advance_ALM_parameters` (default `1.0` for both), `PenaltyContactProcess` multiplies the result by $$10^4$$ afterwards and `MeshTyingProcess` calls it with `compute_penalty: false` to obtain only the scale factor. Test: `tests/cpp_tests/processes/test_alm_variables_calculation_process.cpp` (`ALMVariablesProcess`).

### `AALMAdaptPenaltyValueProcess`

[`aalm_adapt_penalty_value_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/aalm_adapt_penalty_value_process.h) implements the *adapted augmented Lagrangian method* of Bussetta, Marceau and Ponthot (thesis App. D.4.3.1, Algorithm 7), which updates the nodal penalty every iteration from the evolution of the normalized weighted gap $$g_i = \tilde{g}_{n}/A$$ between the current iteration ($$g_i$$, `WEIGHTED_GAP`) and the previous one ($$g_{i-1}$$, `WEIGHTED_GAP` in buffer position 1, stored there by `BaseMortarConvergenceCriteria::PostCriteria`). The thesis algorithm, transcribed from p. 320:

```text
Algorithm 7  Adaptation of normal penalty coefficient [Bussetta, Marceau, Ponthot]
Require: ε_n, g_i and g_{i-1}
 1: procedure ADAPTATION OF NORMAL PENALTY COEFFICIENT
 2:   if g_i × g_{i-1} < 0 then                       -- the gap changed sign
 3:     if g_i × g_{i-1} < 0 then                     -- (sic; the code tests |g_{i-1}| > g_max)
 4:       ε_n = | (ε_n g_{i-1}) / g_i × (|g_i| + g_max) / (g_i − g_{i-1}) |
 5:     else
 6:       ε_n = | ε_n g_{i-1} / (10 g_i) |
 7:   else if g_i > g_max then                        -- penetration above the limit
 8:     if |g_i − g_{i-1}| > max(g_i/10, g_{i-1}/10, 5 g_max) then
 9:       ε_n = 2 ε_n
10:     else if |g_i| = |g_{i-1}| ± 1%  <  10 g_max then
11:       ε_n = ε_n ( sqrt(|g_i|/g_max − 1) + 1 )^2
12:     else if g_i > g_max then
13:       ε_n = 2 ε_n (g_{i-1}/g_i)
14:     else
15:       ε_n = ε_n ( sqrt(|g_i|/g_max − 1) + 1 )
16:   else                                            -- penetration below the limit
17:     ε_n = ε_n
```

In the implementation $$g_{max}$$ = `MAX_GAP_FACTOR` $$\times$$ `NODAL_H` of the node (`MAX_GAP_FACTOR` is set by `ALMContactProcess` from `advance_ALM_parameters["max_gap_factor"]`, default `1.0e-3`), the starting value is the global `INITIAL_PENALTY` in the first iteration of the first step and the nodal `INITIAL_PENALTY` afterwards, absolute values are used wherever the paper compares gaps with the limit (`std::abs(previous_gap) > max_gap`, `std::abs(current_gap) > max_gap`, the comment in the source notes that the abs is "deduced from the paper"), the line-12 test becomes `std::abs(current_gap) > std::abs(previous_gap) * 1.01`, and nodes without `NODAL_AREA` keep their penalty. The result is written to the nodal `INITIAL_PENALTY`, which the conditions and `ActiveSetUtilities` read in preference to the global value. It is enabled with `advance_ALM_parameters["adapt_penalty"] = true` (`ProcessInfo[ADAPT_PENALTY]`) and executed by `BaseMortarConvergenceCriteria::PreCriteria`. Test: `tests/cpp_tests/processes/test_aalm_processes.cpp` (`AALMProcess1`).

### `ComputeDynamicFactorProcess`

[`compute_dynamic_factor_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/compute_dynamic_factor_process.h) stabilizes impact problems. For every node that is `SLAVE` and `ACTIVE` it computes the normalized current and previous gaps $$g_i = \tilde{g}_n/A$$, $$g_{i-1}$$ (writing $$g_i$$ to `NORMAL_GAP` for the post-process) and:

- if `MAX_GAP_THRESHOLD` $$\gt 0$$ and $$g_i \le 0$$, scales the nodal penalty with a logistic function of the penetration, $$\varepsilon_i = \varepsilon \left(1 + f\, \mathrm{MAX\_GAP\_FACTOR}\right)$$ with $$f = 1/(1 + e^{-6 \vert g_i \vert / g_{thr}})$$ (`ComputeLogisticFactor`); otherwise the nodal penalty is reset to the global one;
- if the node just entered contact ($$g_i \lt 0$$ and $$g_{i-1} \gt 0$$), sets `DYNAMIC_FACTOR` $$= \min\left(1, \vert g_i \vert / \vert g_i - g_{i-1} \vert\right)$$, the fraction of the step during which the node was actually in contact. The conditions multiply the contact contribution by this factor (see [Gap history, `Predict()` and the dynamic factor](../Contact_Search/Gap_Computation.html#gap-history-predict-and-the-dynamic-factor)).

It is executed from `BaseMortarConvergenceCriteria::PreCriteria` when the problem is dynamic (`VELOCITY` is a nodal variable) and `compute_dynamic_factor` is set in the solver `contact_settings`, and every step from `ExplicitPenaltyContactProcess.ExecuteInitializeSolutionStep`, which also sets `MAX_GAP_THRESHOLD` (mean `NODAL_H`, or `advance_explicit_parameters["max_gap_threshold"]` when `manual_max_gap_theshold`).

## Model-part bookkeeping processes

### `MasterSlaveProcess`

[`master_slave_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/master_slave_process.h) is used when the `Contact` sub-model-part already exists (restart, or interface conditions supplied by the user) and the `MASTER`/`SLAVE` flags exist only on the nodes. `Execute()` collects every node that is `INTERFACE`, sets each condition whose nodes are all `INTERFACE` to `SLAVE` if all its nodes are `SLAVE` and to `MASTER` otherwise, and adds those nodes and conditions to `Contact`. `SearchBaseProcess.ExecuteInitialize` runs it when `preprocess` is `False`.

### `AssignParentElementConditionsProcess`

[`assign_parent_element_conditions_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/assign_parent_element_conditions_process.h) links every condition of one model part to the element of another model part that owns its face. Defaults (for the `(Model&, Parameters)` constructor):

```json
{
    "conditions_model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
    "elements_model_part_name"   : "PLEASE_SPECIFY_MODEL_PART_NAME",
    "echo_level"                 : 0
}
```

`ExecuteInitialize` builds a hash map from the sorted node ids of every element boundary (`GenerateBoundariesEntities`) to the element id; `ExecuteInitializeSolutionStep` looks up the sorted node ids of every condition and stores the element pointer in the condition variable `PARENT_ELEMENT` (a warning is printed when a face is not found); `Execute` runs both. It is used by `MeshTyingProcess` when `consider_static_condensation` is `true`, so that the mesh-tying conditions can access their parent element. Test: `tests/cpp_tests/processes/test_assign_parent_element_conditions_process.cpp` (`AssignParentElementConditionsProcess1`).

### `FindIntersectedGeometricalObjectsWithOBBContactSearchProcess`

[`find_intersected_geometrical_objects_with_obb_for_contact_search_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/find_intersected_geometrical_objects_with_obb_for_contact_search_process.h) specializes the core `FindIntersectedGeometricalObjectsWithOBBProcess` for the octree broad phase (`type_search: "OctreeWithOBB"`). It works on conditions only, overrides `SetOctreeBoundingBox` and `MarkIfIntersected` (the intersecting masters are flagged `SELECTED`) and supports an asymmetric enlargement of the OBB through `lower_bounding_box_coefficient` / `higher_bounding_box_coefficient` (members `mLowerBBCoefficient = 0.0`, `mHigherBBCoefficient = 1.0`). Defaults:

```json
{
    "intersected_model_part_name"     : "",
    "intersecting_model_part_name"    : "",
    "bounding_box_factor"             : -1.0,
    "debug_obb"                       : false,
    "OBB_intersection_type"           : "SeparatingAxisTheorem",
    "build_from_bounding_box"         : true,
    "lower_bounding_box_coefficient"  : 0.0,
    "higher_bounding_box_coefficient" : 1.0,
    "intersecting_conditions"         : true,
    "intersecting_elements"           : false,
    "intersected_conditions"          : true,
    "intersected_elements"            : false
}
```

`BaseContactSearchProcess::SearchUsingOcTree` constructs it with `MasterSubModelPart<N>` / `SlaveSubModelPart<N>` and the `octree_search_parameters` block (with `bounding_box_factor` multiplied by the maximum `NODAL_H`), then calls `IdentifyNearEntitiesAndCheckEntityForIntersection` for each active slave condition. It is not compatible with `inverted_search`. Test: `tests/cpp_tests/processes/test_search_process.cpp` (`SearchProcessOctree`).

## Error estimation

### `ContactSPRErrorProcess`

[`contact_spr_error_process.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/contact_spr_error_process.h) (author Anna Rehr) extends the MeshingApplication `SPRErrorProcess<TDim>` (superconvergent patch recovery) with the contact constraints: in `CalculatePatch`, for nodes flagged `CONTACT`, the least-squares system of the patch is augmented with penalty terms that enforce the recovered normal stress to match the contact pressure (`AUGMENTED_NORMAL_CONTACT_PRESSURE`, or `CONTACT_PRESSURE` when the former is missing) and the recovered tangential stresses to vanish, using the normal/tangent matrices of `ComputeNormalTangentMatrices` (2D and 3D specializations). Defaults:

```json
{
    "stress_vector_variable"              : "CAUCHY_STRESS_VECTOR",
    "penalty_normal"                      : 1.0e4,
    "penalty_tangential"                  : 1.0e4,
    "echo_level"                          : 0
}
```

It is instantiated as `ContactSPRErrorProcess2D` / `3D` by `ContactRemeshMmgProcess._GenerateErrorProcess` and by the adaptive-remeshing solvers (`compute_error_extra_parameters`), see [Adaptive remeshing](../Examples/Adaptive_Remeshing.html).

## Python processes

The Python processes are the classes the user names in `ProjectParameters.json`. Their complete default JSON blocks and the meaning of every key are in the [Contact process settings reference](../Usage/Contact_Process_Settings_Reference.html); this section documents their inheritance, what each `Execute*` stage does and which C++ objects they drive.

| Module (`python_scripts/`) | Class | Base | Condition name produced | Creates |
|---|---|---|---|---|
| `search_base_process.py` | `SearchBaseProcess` | `KM.Process` | (abstract, `_get_condition_name` returns `""`) | `Contact` model parts, `InterfacePreprocessCondition`, `NormalCheckProcess`, `MasterSlaveProcess`, `CSMA.ContactSearchProcess` per pair |
| `alm_contact_process.py` | `ALMContactProcess` | `SearchBaseProcess` | `ALM[NV]Frictionless[Components][Axisym]MortarContact`, `ALM[NV]Frictional[Axisym]MortarContact` | + `ALMVariablesCalculationProcess`, `ALMFastInit` |
| `penalty_contact_process.py` | `PenaltyContactProcess` | `ALMContactProcess` | `Penalty[NV]Frictionless[Axisym]MortarContact`, `Penalty[NV]Frictional[Axisym]MortarContact` | same as ALM, penalty rescaled |
| `explicit_penalty_contact_process.py` | `ExplicitPenaltyContactProcess` | `PenaltyContactProcess` | same as penalty | + `ComputeDynamicFactorProcess`, `ActiveSetUtilities`, `ContactUtilities` |
| `mpc_contact_process.py` | `MPCContactProcess` | `SearchBaseProcess` | `MPCMortarContact` | `CSMA.MPCContactSearchProcess` per pair, `ALMFastInit` |
| `mesh_tying_process.py` | `MeshTyingProcess` | `SearchBaseProcess` | `MeshTyingMortar` | + `ALMVariablesCalculationProcess` (scale factor only), `AssignParentElementConditionsProcess` |
| `contact_remesh_mmg_process.py` | `ContactRemeshMmgProcess` | `MmgProcess` (MeshingApplication) | – | `ContactSPRErrorProcess2D/3D`, cleans `Contact` before remeshing |
| `basic_mapping_process.py` | `BasicMappingProcess` | `KM.Process` | – | core `SimpleMortarMapperProcess` |
| `replace_properties_process.py` | `ReplacePropertiesProcess` | `KM.Process` | – | core `ReadMaterialsUtility` |

The `[NV]` infix appears when `normal_variation` is `"elemental_derivatives"` (`NODAL_ELEMENTAL_DERIVATIVES`), `[Axisym]` when `alternative_formulations["axisymmetric"]` is `true` (2D only). The C++ search appends `Condition<TDim>D<TNumNodes>N` plus the `final_string` (`_get_final_string`: `"<TNumNodesMaster>N"` when master and slave have a different number of nodes, empty otherwise) to this base name, e.g. `ALMFrictionlessMortarContactCondition3D3N4N`.

### `SearchBaseProcess` life cycle

`SearchBaseProcess` is never used directly; it receives from the derived class a reduced settings block (`search_model_part`, `search_property_ids`, `assume_master_slave`, `model_part_name`, `interval`, `zero_tolerance_factor`, `integration_order`, `consider_tessellation`, `normal_check_proportion`, `search_parameters`) and implements the common machinery. The derived classes customize it through the protected hooks `_get_condition_name`, `_get_final_string`, `_get_problem_name`, `_initialize_process_info`, `_initialize_search_values`, `_initialize_problem_parameters`, `_initialize_search_conditions`, `_set_additional_parameters` and `_create_main_search`.

| Stage | What `SearchBaseProcess` does |
|---|---|
| `__init__` | Validates the settings, stores `main_model_part`, `dimension` (`DOMAIN_SIZE`), `database_step = 0`, `predefined_master_slave`, and the `IntervalUtility` |
| `ExecuteInitialize` | (1) Creates or reuses the `Contact` sub-model-part; if the root model part is `MODIFIED` (remeshing) it first calls `ContactUtilities.CleanContactModelParts`, removes `Contact` and re-creates it. (2) When preprocessing is needed, builds `ContactSub<N>` for every non-empty key of `search_model_part` (`__generate_search_model_part_from_input_list`: `INTERFACE`, `MASTER`, `SLAVE` flags, pair properties `100 + N` or `search_property_ids[N]`, `InterfacePreprocessCondition.GenerateInterfacePart`, `FastTransferBetweenModelPartsProcess`) or detects the skin (`__detect_skin`, `SkinDetectionProcess2D/3D`) when no model part is given. (3) `FindNodalHProcess`. (4) `NormalCheckProcess` unless `IS_RESTARTED`. (5) `adapt_search`: multiplies `search_factor` and `active_check_factor` by `ContactUtilities.CalculateRelativeSizeMesh`. (6) `_initialize_process_info` (`ZERO_TOLERANCE_FACTOR`, `ACTIVE_CHECK_FACTOR`). (7) `MasterSlaveProcess` when the conditions already existed. (8) `INTEGRATION_ORDER_CONTACT` and `CONSIDER_TESSELLATION` on every property of `Contact`. (9) `_initialize_search_values` (`DISTANCE_THRESHOLD = 1.0e24`, `ACTIVE_CHECK_FACTOR` on the properties). (10) `_initialize_problem_parameters`. (11) `_create_main_search` per pair. (12) `_initialize_search_conditions`. (13) `ExecuteInitialize` of every C++ search |
| `ExecuteBeforeSolutionLoop` | nothing |
| `ExecuteInitializeSolutionStep` | If `_compute_search()` (inside the interval and `STEP == 1` or `database_step >= database_step_update`): resets `MARKER` on the nodes, calls `ExecuteInitializeSolutionStep` of every C++ search (optionally `_debug_output` GiD dump per pair), resets `MARKER` again and, in `debug_mode`, prints the integrated contact area (`__get_integration_area`, `ExactMortarIntegrationUtility*`). Otherwise sets `ACTIVE = false` on all nodes and conditions of `Contact` |
| `ExecuteFinalizeSolutionStep` | Inside the interval, if the model part is not `MODIFIED` and the search was run this step (same counter test), calls `ExecuteFinalizeSolutionStep` of every C++ search (`ClearMortarConditions`) and resets `database_step` |
| `ExecuteBeforeOutputStep`, `ExecuteAfterOutputStep`, `ExecuteFinalize` | nothing |

`_create_search_parameters` builds the C++ `Parameters` from `search_parameters`: it copies `simple_search`, `type_search`, `check_gap`, `bucket_size`, `search_factor`, `dynamic_search`, `static_check_movement`, `consider_gap_threshold`, `normal_orientation_threshold` and `debug_mode`, maps `max_number_results` to `allocation_size`, fills `condition_name`, `final_string`, `predefined_master_slave` (`false` when `assume_master_slave[N]` is empty) and `id_name = N`, and lets the derived class add keys with `_set_additional_parameters` (the contact processes add `pure_slip`). Note that `active_check_factor`, `database_step_update`, `adapt_search`, `predict_correct_lagrange_multiplier` and `octree_search_parameters` are consumed on the Python side or passed through the `ProcessInfo`, not copied into the C++ parameters.

### `ALMContactProcess`

`__init__` validates the full ALM block, transfers the shared keys to the base (`contact_model_part` → `search_model_part`, `contact_property_ids` → `search_property_ids`), resolves `normal_variation` into a `CSMA.NormalDerivativesComputation` value, stores `frictional_law`, and derives `is_frictional`, `pure_slip` and `slip_step_reset_frequency` from `contact_type` (`"Frictionless"`, `"FrictionlessComponents"`, `"Frictional"`, with optional `PureSlip` / `WithNormalUpdate` suffixes; frictional problems get `WithNormalUpdate` unless `not_normal_update_frictional`). 3D axisymmetric input raises an error.

| Stage | Additional work of `ALMContactProcess` |
|---|---|
| `ExecuteInitialize` | Frictional case: `pure_slip` is `True` for `PureSlip` types, or when all `friction_coefficients` are zero and `auxiliary_methods_solvers.AuxiliaryPureSlipCheck` confirms it. Zeroes `AUGMENTED_NORMAL_CONTACT_PRESSURE` (and `AUGMENTED_TANGENT_CONTACT_PRESSURE`) on all nodes, then calls the base |
| `_initialize_process_info` | `CONSIDER_NORMAL_VARIATION`, `ACTIVE_SET_CONVERGED = True`, `ADAPT_PENALTY`, `MAX_GAP_FACTOR`, `OPERATOR_THRESHOLD` |
| `_initialize_search_values` | `CONTACT` flag on the main and `Contact` model parts, `SLIP` flag = `is_frictional`; frictional: `TANGENT_FACTOR`, `SLIP_AUGMENTATION_COEFFICIENT`, `SLIP_THRESHOLD` |
| `_initialize_problem_parameters` | Unless `manual_ALM`: `ALMVariablesCalculationProcess(Contact, NODAL_H, {stiffness_factor, penalty_scale_factor})`, then `SCALE_FACTOR = 1` if `use_scale_factor` is `false`; with `manual_ALM`: `INITIAL_PENALTY = penalty`, `SCALE_FACTOR = scale_factor`. Both are clamped to at least `1.0` and printed |
| `_initialize_search_conditions` | Frictional: writes `friction_coefficients[N]` into `FRICTION_COEFFICIENT` of the properties of `ContactSub<N>` (warns if already present). Then `ALMFastInit(Contact).Execute()` |
| `_set_additional_parameters` | Adds `pure_slip` to the C++ search parameters |
| `ExecuteInitializeSolutionStep` | Base search, then `_reset_slip_flag`: every `slip_step_reset_frequency` steps the nodal `SLIP` flags are cleared (`> 0`); `0` never resets; `< 0` recomputes the tangent from `WEIGHTED_SLIP` with `MortarUtilities.ComputeNodesTangentModelPart` |
| `ExecuteFinalizeSolutionStep` | Base; in `debug_mode` prints the total applied load (`LINE_LOAD`/`SURFACE_LOAD`), total `REACTION` and total contact force $$\sum A_i \bar{\lambda}_{n,i}$$ |
| `ExecuteBeforeOutputStep` | With `clear_inactive_for_post` (default `true`) zeroes `AUGMENTED_NORMAL_CONTACT_PRESSURE` and `AUGMENTED_TANGENT_CONTACT_PRESSURE` on the nodes that are not `ACTIVE`, so that the post-process shows pressure only on the active set |

### `PenaltyContactProcess`

Same settings as `ALMContactProcess` except `tangent_factor: 1.0e-3` and, in `advance_ALM_parameters`, `penalty: 1.0e16` and `max_gap_factor: 5.0e-4`. It overrides `_get_condition_name` (the `Penalty*MortarContact` names) and `_initialize_problem_parameters`: unless `manual_ALM`, it recomputes `NODAL_H`, runs `ALMVariablesCalculationProcess` and multiplies the resulting `INITIAL_PENALTY` by $$10^4$$ ("the process is designed for the ALM formulation"); the minimum penalty is `1.0e16` instead of `1.0`. No scale factor is used by the penalty conditions.

### `ExplicitPenaltyContactProcess`

Designed for `contact_explicit_dynamic_solver`. Defaults differ from the penalty block in `tangent_factor: 1.0e-4`, `search_parameters.type_search: "octree_with_obb"` and `advance_ALM_parameters.max_gap_factor: 1.0e-3`.

| Stage | Additional work |
|---|---|
| `ExecuteInitialize` | Base; then `NL_ITERATION_NUMBER = 1` (the active-set utilities test it), `MAX_GAP_THRESHOLD` = `max_gap_threshold` (when `manual_max_gap_theshold`) or the mean `NODAL_H` (`ContactUtilities.CalculateMeanNodalH`), and creates `ComputeDynamicFactorProcess` |
| `ExecuteInitializeSolutionStep` | Base search; `ContactUtilities.CheckActivity(main, False)` sets the `CONTACT` flag of the main model part (the explicit solver reduces the time step by `delta_time_factor_for_contact` when it is set). If active and inside the interval: zeroes `NODAL_AREA` and `WEIGHTED_GAP`, `ContactUtilities.ComputeExplicitContributionConditions`, `ActiveSetUtilities.ComputePenaltyFrictionlessActiveSet` or `ComputePenaltyFrictionalActiveSet`, `ContactUtilities.ActivateConditionWithActiveNodes` and `ComputeDynamicFactorProcess.Execute()` |
| `_compute_search` | Same counter logic as the base (a comment in the source anticipates an explicit-specific criterion) |
| `_initialize_problem_parameters` | As the penalty process (its own copy, to avoid computing the values twice) |

Since there is no Newton loop in explicit dynamics, this process performs itself what the convergence criteria do in the implicit solvers (explicit contribution, active set, dynamic factor).

### `MPCContactProcess`

Constraint-based contact (see [Frictional laws and MPC constraint](Frictional_Laws_And_MPC_Constraint.html)). Its defaults share the ALM structure but with `tangent_factor: 1.0e-1`, `zero_tolerance_factor: 1.0e2`, `reaction_check_stiffness_factor: 1.0e-10`, `update_condition_relation_step: false` and without `advance_ALM_parameters` / `alternative_formulations` / slip keys.

| Stage | Additional work |
|---|---|
| `ExecuteInitialize` | Same `pure_slip` detection as the ALM process, then the base |
| `_initialize_process_info` | `ACTIVE_SET_CONVERGED = True`, `REACTION_CHECK_STIFFNESS_FACTOR` |
| `_initialize_search_values` | `CONTACT` and `SLIP` flags; frictional: `TANGENT_FACTOR` |
| `_initialize_search_conditions` | Friction coefficients on the pair properties, `ALMFastInit` |
| `_create_main_search` | `CSMA.MPCContactSearchProcess(main_model_part, search_parameters)` (no pair properties argument) |
| `ExecuteInitializeSolutionStep` | Base search; frictional: sets `SLIP` on all conditions |
| `ExecuteFinalizeSolutionStep` | Base; unless `update_condition_relation_step`, sets `BLOCKED` on all conditions so that the constraint relation is frozen for the step |

Note: `mpc_contact_process.py` defines `ExecuteFinalizeSolutionStep` twice; the second definition (the `BLOCKED` logic) is the one in effect and the first one (a `debug_mode` load/reaction balance identical to the ALM one) is unreachable.

### `MeshTyingProcess`

Mortar mesh tying ([Mesh tying](../Theory/Mesh_Tying.html)). Its settings rename the pair keys (`mesh_tying_model_part`, `mesh_tying_property_ids`), add `variable_name` (`"DISPLACEMENT"`; a scalar variable gives the `Scalar` tying, a vector one the `Components` tying), `consider_static_condensation`, `scale_factor_parameters` (`manual_scale_factor`, `stiffness_factor`, `scale_factor`), default `consider_tessellation: true` and `database_step_update: 999999999` (the pairing is computed once).

| Stage | Additional work |
|---|---|
| `ExecuteInitialize` | Base; then `scale_factor` from `scale_factor_parameters` or `ALMVariablesCalculationProcess(Contact, NODAL_H, {compute_penalty: false, stiffness_factor})`; with `consider_static_condensation`, creates `AssignParentElementConditionsProcess(ComputingContact, main_model_part)` and runs its `ExecuteInitialize` |
| `ExecuteInitializeSolutionStep` | Base; with static condensation, `AssignParentElementConditionsProcess.ExecuteInitializeSolutionStep` (re-links the new paired conditions to their parent elements) |
| `_get_condition_name` | `"MeshTyingMortar"` |
| `_initialize_search_conditions` | `TYING_VARIABLE = variable_name` on the properties; zeroes `WEIGHTED_SCALAR_RESIDUAL` or `WEIGHTED_VECTOR_RESIDUAL` |

### `ContactRemeshMmgProcess`

Derives from `MmgProcess` of the MeshingApplication and adapts it to contact: the Hessian metric is computed on `VON_MISES_STRESS`, `AUGMENTED_NORMAL_CONTACT_PRESSURE` and optionally `STRAIN_ENERGY` (`consider_strain_energy`), with an automatic normalization factor $$20/(\nu^2 E)$$ (`automatic_normalization_factor`) read from the first properties with `YOUNG_MODULUS` and `POISSON_RATIO`; `_AuxiliaryCallsBeforeRemesh` calls `ContactUtilities.CleanContactModelParts`, removes all entities flagged `TO_ERASE` (including master–slave constraints) and re-creates the empty `Contact` sub-model-part so that the next `SearchBaseProcess.ExecuteInitialize` (triggered by the `MODIFIED` flag) rebuilds the interface; `_AuxiliaryCallsAfterRemesh` transfers the new entities to the parent model part; `_GenerateErrorProcess` returns `ContactSPRErrorProcess2D/3D` configured from `error_strategy_parameters["compute_error_extra_parameters"]` (`stress_vector_variable`, `penalty_normal`, `penalty_tangential`). It requires the MeshingApplication compiled with MMG; see [Adaptive remeshing](../Examples/Adaptive_Remeshing.html).

### `BasicMappingProcess`

A thin `KM.Process` around the core `SimpleMortarMapperProcess`: the settings are the mapper settings plus `origin_model_part_name`, `destination_model_part_name`, `interval` (default `[0.0, 1e30]`) and an optional `linear_solver_settings` block (a solver is built with `linear_solver_factory` when `solver_type` is given). `ExecuteInitializeSolutionStep` runs the mapper when `TIME` is inside the interval. It is useful to transfer results between non-matching meshes in the same input file that uses the contact processes.

### `ReplacePropertiesProcess`

Reloads the materials of `model_part_name` from `materials_filename` with `ReadMaterialsUtility` at every step whose `TIME` lies inside `interval` and, when `reinitialize_entities` is `true`, calls `Initialize` on every element and condition of the model part. The `"End"` string is accepted as the second interval value.

Note: the default-settings string of `replace_properties_process.py` contains `"reinitialize_entities" : false.` — a period instead of a comma before `"interval"` — so `KM.Parameters(...)` raises a JSON parse error and the process cannot be constructed in its present form. It is not used by any test or example of the application.

## Call sequence

The following text diagrams show who calls whom. `SBP` stands for the Python process (`SearchBaseProcess` and its derived class), `CSP` for the C++ search (`ContactSearchWrapperProcess` → `Simple`/`AdvancedContactSearchProcess`, or the MPC variant).

**Initialization** (`AnalysisStage.Initialize` → `ExecuteInitialize` of every process in `contact_process_list`):

```text
SBP.ExecuteInitialize
├── [ALM/Penalty/MPC] pure-slip detection, zero AUGMENTED_*_CONTACT_PRESSURE
├── create/reuse "Contact"   (MODIFIED → ContactUtilities.CleanContactModelParts, rebuild)
├── for every pair N: __generate_search_model_part_from_input_list
│   ├── flags INTERFACE / MASTER / SLAVE (_assign_master_flags, _assign_slave_flags)
│   ├── InterfacePreprocessCondition.GenerateInterfacePart   (nodes → face conditions)
│   └── FastTransferBetweenModelPartsProcess → "Contact/ContactSubN"
├── FindNodalHProcess
├── NormalCheckProcess.Execute                               (unless IS_RESTARTED)
├── [adapt_search] ContactUtilities.CalculateRelativeSizeMesh
├── _initialize_process_info      → ZERO_TOLERANCE_FACTOR, ACTIVE_CHECK_FACTOR, [ALM] CONSIDER_NORMAL_VARIATION, ADAPT_PENALTY, MAX_GAP_FACTOR, OPERATOR_THRESHOLD
├── MasterSlaveProcess.Execute                               (only if "Contact" pre-existed)
├── _initialize_search_values     → DISTANCE_THRESHOLD, CONTACT / SLIP flags, TANGENT_FACTOR, ...
├── _initialize_problem_parameters→ ALMVariablesCalculationProcess.Execute  → INITIAL_PENALTY, SCALE_FACTOR
├── for every pair N: _create_main_search → CSP = CSMA.ContactSearchProcess(main, params, pair properties)
│   └── BaseContactSearchProcess ctor: "ComputingContact/ComputingContactSubN", prototype condition, TypeSolution
├── _initialize_search_conditions → FRICTION_COEFFICIENT on pair properties, ALMFastInit.Execute
└── for every pair N: CSP.ExecuteInitialize
    ├── CheckContactModelParts
    ├── CreatePointListMortar
    └── InitializeMortarConditions
```

**Beginning of a time step** (`AnalysisStage.InitializeSolutionStep` → `ExecuteInitializeSolutionStep`; in the implicit solvers the same call is repeated by `ResidualBasedNewtonRaphsonContactStrategy::AdaptativeStep` through `ProcessFactoryUtility` when a step is split):

```text
SBP.ExecuteInitializeSolutionStep
├── _compute_search()  →  false: ACTIVE = false on Contact nodes and conditions; return
├── reset MARKER on Contact nodes
├── for every pair N: CSP.ExecuteInitializeSolutionStep
│   ├── ClearMortarConditions
│   │   ├── ResetContactOperators          (remove inactive paired conditions [+ MPC constraints])
│   │   └── zero the multiplier of inactive nodes (by TypeSolution)
│   └── UpdateMortarConditions
│       ├── UpdatePointListMortar          ([dynamic_search] ContactUtilities.ComputeStepJump)
│       ├── unit normals of ContactSubN; [not predefined] SelfContactUtilities.NotPredefinedMasterSlave
│       ├── SearchUsingKDTree  |  SearchUsingOcTree → FindIntersectedGeometricalObjectsWithOBBContactSearchProcess
│       │   └── per candidate: OBB test, CheckGeometricalObject → INDEX_MAP  (or AddPotentialPairing / AddPairing)
│       ├── [MappingCheck] CheckPairing
│       │   ├── ComputeMappedGap → NormalGapProcess(MasterSubN, SlaveSubN).Execute → SimpleMortarMapperProcess
│       │   ├── CreateAuxiliaryConditions → AddPairing (clone prototype, ComputingContactSubN)
│       │   └── ComputeWeightedReaction  → ContactUtilities::ComputeExplicitContributionConditions
│       └── ComputeActiveInactiveNodes → SetActiveNode / SetInactiveNode  (Simple or Advanced)
├── reset MARKER; [debug_mode] _debug_output, __get_integration_area
├── [ALM] _reset_slip_flag
├── [MPC] SLIP on conditions
└── [Explicit] CheckActivity → CONTACT flag; ComputeExplicitContributionConditions;
              ActiveSetUtilities.ComputePenalty*ActiveSet; ActivateConditionWithActiveNodes; ComputeDynamicFactorProcess.Execute
```

**Inside the Newton loop** (implicit solvers). The processes are not called directly; the mortar convergence criteria created by `contact_convergence_criteria_factory.py` call them (see [Strategies and convergence criteria](Strategies_And_Convergence_Criteria.html)):

```text
MortarAndConvergenceCriteria
└── <Formulation>MortarConvergenceCriteria : BaseMortarConvergenceCriteria
    ├── PreCriteria  (before the residual is evaluated)
    │   ├── [normal variation] ComputeNodesMeanNormalModelPartWithPairedNormal
    │   ├── [frictional] MortarUtilities::ComputeNodesTangentModelPart
    │   ├── [ADAPT_PENALTY or dynamic] ResetWeightedGap + ContactUtilities::ComputeExplicitContributionConditions
    │   ├── [dynamic and COMPUTE_DYNAMIC_FACTOR] ComputeDynamicFactorProcess(Contact).Execute
    │   └── [ADAPT_PENALTY] AALMAdaptPenaltyValueProcess(Contact).Execute
    └── PostCriteria (after the linear solve)
        ├── WEIGHTED_GAP → buffer 1; ResetWeightedGap; ComputeExplicitContributionConditions
        └── ActiveSetUtilities::Compute<Formulation>ActiveSet  → ACTIVE / SLIP flags, ACTIVE_SET_CONVERGED
```

**End of the step**: `SBP.ExecuteFinalizeSolutionStep` → `CSP.ExecuteFinalizeSolutionStep` → `ClearMortarConditions` (the inactive pairs are deleted; the active ones survive to the next step); `[MPC]` `BLOCKED` on conditions. `SBP.ExecuteBeforeOutputStep` → `[ALM]` clear pressures on inactive nodes.

**`ProcessFactoryUtility`.** The contact solvers wrap the Python process lists in `CSMA.ProcessFactoryUtility` (`AddProcessesList`, `AddPostProcess`) and hand them to `ResidualBasedNewtonRaphsonContactStrategy` / `LineSearchContactStrategy`, whose `AdaptativeStep` re-executes `ExecuteInitializeSolutionStep`, `ExecuteFinalizeSolutionStep`, `ExecuteBeforeOutputStep`, `PrintOutput` and `ExecuteAfterOutputStep` on every sub-step when `adaptative_strategy` is enabled. See [Utilities](Utilities.html#processfactoryutility).

## Tests

| Test | Location | Processes covered |
|---|---|---|
| `SearchProcessKDTree`, `SearchProcessKDTreeWithOBB`, `SearchProcessOctree` | `tests/cpp_tests/processes/test_search_process.cpp` | `ContactSearchWrapperProcess` (→ `AdvancedContactSearchProcess`), `FindIntersectedGeometricalObjectsWithOBBContactSearchProcess` |
| `WeightedGap1` … `WeightedGap9` (11 cases) | `tests/cpp_tests/processes/test_weighted_gap.cpp` | Pairing and explicit contribution of the search output (`AddExplicitContribution`) |
| `AALMProcess1` | `tests/cpp_tests/processes/test_aalm_processes.cpp` | `AALMAdaptPenaltyValueProcess` |
| `ALMVariablesProcess` | `tests/cpp_tests/processes/test_alm_variables_calculation_process.cpp` | `ALMVariablesCalculationProcess` |
| `AssignParentElementConditionsProcess1` | `tests/cpp_tests/processes/test_assign_parent_element_conditions_process.cpp` | `AssignParentElementConditionsProcess` |
| `test_dynamic_search_triangle`, `test_dynamic_search_quad` | `tests/test_dynamic_search.py` | `InterfacePreprocessCondition`, `ALMFastInit`, `ContactSearchProcess` with `dynamic_search` |
| `test_check_normals`, `test_check_normals_quads`, `test_check_normals_s_shape` | `tests/test_check_normals_process.py` | `NormalCheckProcess` |
| `test_process_factory`, `test_processes_list_factory` | `tests/test_process_factory.py` | `ProcessFactoryUtility` |
| Small/Nightly/Validation suites | `tests/SmallTests.py`, `tests/NightlyTests.py`, `tests/ValidationTests.py` | End-to-end runs of `ALMContactProcess`, `PenaltyContactProcess`, `ExplicitPenaltyContactProcess`, `MPCContactProcess`, `MeshTyingProcess` (see [Test suite reference](../Validation/Test_Suite_Reference.html)) |

## Notes and limitations

- `NormalGapProcess` and `AALMAdaptPenaltyValueProcess` are building blocks of the search and of the criteria respectively; running them from `ProjectParameters.json` is not supported (the former needs `AUXILIAR_COORDINATES`, `DISTANCE_THRESHOLD` and the master/slave sub-model-parts prepared by the search; the latter is not exposed to Python).
- `BaseContactSearchProcess` requires a `Contact` sub-model-part and, with `id_name`, the `ContactSub<id_name>` sub-model-part; the constructor raises otherwise. `SearchBaseProcess` always creates them, but direct use from Python (as in `tests/test_dynamic_search.py`) must reproduce that structure.
- The key `manual_max_gap_theshold` of `advance_explicit_parameters` is spelled that way in the source (missing `r`); the same spelling must be used in the input file.
- `MPCContactProcess` has a duplicated `ExecuteFinalizeSolutionStep` (see above) and `ReplacePropertiesProcess` has a malformed default JSON; both are documented as notes and left untouched here.
- `SearchTreeType::Kdop` is declared but not implemented; `ContactSearchWrapperProcess` accepts the `"KDOP"` string and the base constructor raises `KDOP contact search: Not yet implemented`.
- Increasing `database_step_update` above `1` disables contact (all nodes and conditions `ACTIVE = false`) in the steps where the search is skipped; it is intended for mesh tying, where the pairing does not change.
