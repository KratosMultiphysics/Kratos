# Processes

C++ processes of the application. The central ones implement the **contact search** (creation of the slave–master pair conditions in the `ComputingContact` sub-model-part); the others initialise or update the quantities the mortar conditions and the active-set criteria rely on.

![Contact search pipeline](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/images/csma_search_pipeline.svg)

## Search processes

| Header | Class | Python name | Role |
|---|---|---|---|
| `base_contact_search_process.h` | `BaseContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>` | – | Core of the search: point list of the destination conditions, KD-tree (`in_radius`, `in_box`, with/without OBB) or octree search, oriented-bounding-box filtering, gap check (`no_check` / `direct_check` / `mapping_check`), creation and clean-up of the `PairedCondition`s. Enums `SearchTreeType`, `CheckGap`, `TypeSolution`; local flags `INVERTED_SEARCH`, `CREATE_AUXILIAR_CONDITIONS`, `MULTIPLE_SEARCHS`, `PREDEFINE_MASTER_SLAVE`, `PURE_SLIP`. |
| `simple_contact_search_process.h` | `SimpleContactSearchProcess<…>` | `SimpleContactSearchProcess<geometry>` | Activation of the slave nodes by the gap alone (`active_check_factor`). |
| `advanced_contact_search_process.h` | `AdvancedContactSearchProcess<…>` | `AdvancedContactSearchProcess<geometry>` | Activation with prediction/correction of the Lagrange multipliers from a gap–pressure linear regression (`predict_correct_lagrange_multiplier`); default of the Python processes. |
| `contact_search_wrapper_process.h` | `ContactSearchWrapperProcess` | `ContactSearchProcess` | Dimension/geometry-agnostic wrapper that instantiates the right template (`2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N`) from `DOMAIN_SIZE` and the condition geometries. |
| `mpc_contact_search_process.h`, `mpc_contact_search_wrapper_process.h` | `MPCContactSearchProcess<…>`, `MPCContactSearchWrapperProcess` | `MPCContactSearchProcess<geometry>`, `MPCContactSearchProcess` | Same search but creating `MPCMortarContactCondition`s linked to `ContactMasterSlaveConstraint`s. |
| `find_intersected_geometrical_objects_with_obb_for_contact_search_process.h` | `FindIntersectedGeometricalObjectsWithOBBContactSearchProcess` | same | Octree broad phase with oriented bounding boxes adapted to conditions (`OBB_intersection_type`, `bounding_box_factor`, asymmetric coefficients). |

Default search parameters (JSON of `BaseContactSearchProcess`, exposed to the user through `search_parameters` of the Python processes):

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

## Auxiliary processes

| Header | Class | Python name | Role |
|---|---|---|---|
| `normal_gap_process.h` | `NormalGapProcess<…>` | `NormalGapProcess<geometry>` | Nodal `NORMAL_GAP` between master and slave through the mortar mapper (the `mapping_check` gap check). |
| `normal_check_process.h` | `NormalCheckProcess` | same | Detects and repairs inverted / inconsistent condition normals (`length_proportion`, `check_threshold`). |
| `master_slave_process.h` | `MasterSlaveProcess` | same | Assigns the `MASTER` / `SLAVE` flags to conditions and nodes. |
| `alm_fast_init_process.h` | `ALMFastInit` | same | Fast initialisation of the ALM data: `SLIP` flag, nodal `INITIAL_PENALTY`, `WEIGHTED_GAP`, `WEIGHTED_SLIP`, `DYNAMIC_FACTOR`, augmented pressures, condition `NORMAL`. |
| `alm_variables_calculation_process.h` | `ALMVariablesCalculationProcess` | same | Automatic penalty $\varepsilon$ and scale factor $k$ from `YOUNG_MODULUS` and `NODAL_H` (thesis eq. 4.11, $\varepsilon = k \approx$ `stiffness_factor` $E_{mean}/h_{mean}$). |
| `aalm_adapt_penalty_value_process.h` | `AALMAdaptPenaltyValueProcess` | – (used by the criteria) | Adapted augmented Lagrangian: updates the nodal penalty from the gap evolution (Bussetta–Marceau–Ponthot, thesis Algorithm 7), enabled with `adapt_penalty`. |
| `compute_dynamic_factor_process.h` | `ComputeDynamicFactorProcess` | same | `DYNAMIC_FACTOR` from the ratio of current and previous weighted gap (dynamic contact). |
| `assign_parent_element_conditions_process.h` | `AssignParentElementConditionsProcess` | same | Links each condition to its parent element (`NEIGHBOUR_ELEMENTS`); used by the static condensation of mesh tying. |
| `contact_spr_error_process.h` | `ContactSPRErrorProcess<TDim>` | `ContactSPRErrorProcess2D`, `ContactSPRErrorProcess3D` | Superconvergent-patch-recovery error estimator with contact penalty terms, for adaptive remeshing (`ContactErrorMeshCriteria`). |

The Python counterparts that orchestrate these processes (`SearchBaseProcess`, `ALMContactProcess`, `PenaltyContactProcess`, `ExplicitPenaltyContactProcess`, `MPCContactProcess`, `MeshTyingProcess`) live in [`../python_scripts/`](../python_scripts/README.md).

## Full documentation

- [Processes (implementation reference)](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Processes.html) · [source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Processes.md)
- Search algorithms: [search pipeline and bounding volumes](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/Search_Pipeline_And_Bounding_Volumes.md), [gap computation](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/Gap_Computation.md), [self-contact](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/Self_Contact.md)
- User-facing JSON: [contact process settings](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Contact_Process_Settings_Reference.md)
