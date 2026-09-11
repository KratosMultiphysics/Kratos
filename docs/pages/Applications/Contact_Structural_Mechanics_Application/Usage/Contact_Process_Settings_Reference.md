---
title: Contact Process Settings Reference
keywords: alm_contact_process, penalty_contact_process, mpc_contact_process, mesh_tying_process, search_parameters, advance_ALM_parameters, contact_model_part, assume_master_slave
tags: [usage, settings, processes, JSON, reference, search_parameters]
sidebar: contact_structural_mechanics_application
summary: Key-by-key reference of the Python contact processes (ALM, penalty, explicit penalty, MPC and mesh tying) placed in processes.contact_process_list, with the theory quantity and the C++ object behind every parameter.
---

> **Sources.** [`python_scripts/alm_contact_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/python_scripts/alm_contact_process.py), [`penalty_contact_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/python_scripts/penalty_contact_process.py), [`explicit_penalty_contact_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/python_scripts/explicit_penalty_contact_process.py), [`mpc_contact_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/python_scripts/mpc_contact_process.py), [`mesh_tying_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/python_scripts/mesh_tying_process.py), all built on [`search_base_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/python_scripts/search_base_process.py); thesis §4.3.3.3 (calibration of $$k$$, $$\varepsilon$$), §4.3.4 (frictional parameters), §4.4 (search).

## Where the process goes

Every contact simulation has exactly one process per interface family in `processes.contact_process_list`. The process builds the contact sub-model-parts, creates the search utilities, sets the parameters that the conditions read from the `ProcessInfo` and the properties, and writes the matching `mortar_type` into `solver_settings.contact_settings` (see [Solver settings](Solver_Settings_Reference.html)).

```json
"processes" : {
    "contact_process_list" : [{
        "python_module" : "alm_contact_process",
        "kratos_module" : "KratosMultiphysics.ContactStructuralMechanicsApplication",
        "process_name"  : "ALMContactProcess",
        "Parameters"    : {
            "model_part_name"     : "Structure",
            "contact_model_part"  : { "0" : ["Contact_Part_1", "Contact_Part_2"] },
            "assume_master_slave" : { "0" : ["Parts_Parts_Auto2"] },
            "contact_type"        : "Frictionless"
        }
    }]
}
```

| `python_module` | `process_name` | Formulation | Sets `mortar_type` |
|---|---|---|---|
| `alm_contact_process` | `ALMContactProcess` | augmented Lagrangian mortar contact | `ALMContactFrictionless`, `ALMContactFrictionlessComponents`, `ALMContactFrictional[PureSlip]` |
| `penalty_contact_process` | `PenaltyContactProcess` | penalty mortar contact | `PenaltyContactFrictionless`, `PenaltyContactFrictional[PureSlip]` |
| `explicit_penalty_contact_process` | `ExplicitPenaltyContactProcess` | penalty contact for the explicit solver | idem |
| `mpc_contact_process` | `MPCContactProcess` | multipoint-constraint contact | – (uses `mpc_contact_settings`) |
| `mesh_tying_process` | `MeshTyingProcess` | mortar mesh tying | `ScalarMeshTying`, `ComponentsMeshTying` |

## Defining the interfaces: the pair dictionaries

All processes share the same way of declaring interfaces (inherited from `SearchBaseProcess`): dictionaries with the keys `"0"` … `"9"`, one entry per interface (pair of potentially contacting surfaces). Only non-empty entries are used, and each entry creates its own `ContactSearchProcess`, sub-model-parts `ContactSub<key>`, `MasterSubModelPart<key>`, `SlaveSubModelPart<key>` and computing conditions in `ComputingContact`.

| Key | Meaning |
|---|---|
| `contact_model_part` (ALM/penalty/MPC) / `mesh_tying_model_part` (tying) | List of sub-model-part names whose **conditions** form the interface `<key>`. Usually two surfaces (slave and master), but a single self-contacting surface is also valid. If every entry is empty, the skin of the whole model part is detected automatically (`__detect_skin`). |
| `assume_master_slave` | For interface `<key>`, the sub-model-parts that play the **master** role; everything else in the interface is slave. An empty list activates the automatic (self-contact) master/slave assignment, see [Self contact](../Contact_Search/Self_Contact.html). The slave side is where the mortar integration takes place, so the finer / more curved surface is usually the better slave. |
| `contact_property_ids` / `mesh_tying_property_ids` / `search_property_ids` | Property id used for the conditions created for interface `<key>` (`0` = create a new property copied from the elements). The condition parameters (`INTEGRATION_ORDER_CONTACT`, `CONSIDER_TESSELLATION`, `ACTIVE_CHECK_FACTOR`, `FRICTION_COEFFICIENT`) are written into it. |
| `friction_coefficients` (ALM/penalty/MPC) | Coulomb coefficient $$\mu$$ of interface `<key>`, written as `FRICTION_COEFFICIENT` into the pair property (a value already present in the property is kept, with a warning). |

Multi-interface example (two independent contacts, the second frictional-ready):

```json
"contact_model_part"   : { "0" : ["Contact_Punch", "Contact_Blank_Top"], "1" : ["Contact_Blank_Bottom", "Contact_Die"] },
"assume_master_slave"  : { "0" : ["Contact_Punch"],                      "1" : ["Contact_Die"] },
"friction_coefficients": { "0" : 0.1,                                     "1" : 0.2 }
```

## `alm_contact_process` — full defaults

```json
{
    "help"                          : "This class is used in order to compute the contact using a mortar ALM formulation. This class constructs the model parts containing the contact conditions and initializes parameters and variables related with the contact. The class creates search utilities to be used to create the contact pairs",
    "model_part_name"               : "Structure",
    "contact_model_part"            : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
    "assume_master_slave"           : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
    "contact_property_ids"          : {"0": 0,"1": 0,"2": 0,"3": 0,"4": 0,"5": 0,"6": 0,"7": 0,"8": 0,"9": 0},
    "friction_coefficients"         : {"0": 0.0,"1": 0.0,"2": 0.0,"3": 0.0,"4": 0.0,"5": 0.0,"6": 0.0,"7": 0.0,"8": 0.0,"9": 0.0},
    "contact_type"                  : "Frictionless",
    "not_normal_update_frictional"  : false,
    "interval"                      : [0.0,"End"],
    "normal_variation"              : "no_derivatives_computation",
    "frictional_law"                : "Coulomb",
    "tangent_factor"                : 2.5e-2,
    "operator_threshold"            : 1.0e-3,
    "slip_augmentation_coefficient" : 0.0,
    "slip_threshold"                : 2.0e-2,
    "zero_tolerance_factor"         : 1.0,
    "integration_order"             : 2,
    "consider_tessellation"         : false,
    "normal_check_proportion"       : 0.1,
    "clear_inactive_for_post"       : true,
    "slip_step_reset_frequency"     : 1,
    "search_parameters"             : {
        "type_search"                         : "in_radius_with_obb",
        "simple_search"                       : false,
        "adapt_search"                        : false,
        "search_factor"                       : 3.5,
        "active_check_factor"                 : 0.01,
        "max_number_results"                  : 1000,
        "bucket_size"                         : 4,
        "dynamic_search"                      : false,
        "static_check_movement"               : false,
        "database_step_update"                : 1,
        "normal_orientation_threshold"        : 1.0e-1,
        "consider_gap_threshold"              : false,
        "debug_mode"                          : false,
        "predict_correct_lagrange_multiplier" : false,
        "check_gap"                           : "check_mapping",
        "octree_search_parameters" : {
            "bounding_box_factor"             : 0.1,
            "debug_obb"                       : false,
            "OBB_intersection_type"           : "SeparatingAxisTheorem",
            "build_from_bounding_box"         : true,
            "lower_bounding_box_coefficient"  : 0.0,
            "higher_bounding_box_coefficient" : 1.0
        }
    },
    "advance_explicit_parameters"  : {
        "manual_max_gap_theshold"  : false,
        "automatic_gap_factor"     : 1.0e-1,
        "max_gap_threshold"        : 5.0e-2,
        "max_gap_factor"           : 1.0e2,
        "logistic_exponent_factor" : 6.0
    },
    "advance_ALM_parameters" : {
        "manual_ALM"                  : false,
        "stiffness_factor"            : 1.0,
        "penalty_scale_factor"        : 1.0,
        "use_scale_factor"            : true,
        "penalty"                     : 1.0e-12,
        "scale_factor"                : 1.0e0,
        "adapt_penalty"               : false,
        "max_gap_factor"              : 1.0e-3
    },
    "alternative_formulations" : {
        "axisymmetric"                : false
    }
}
```

### Formulation switches

| Key | Values / default | Meaning |
|---|---|---|
| `contact_type` | `Frictionless` (default), `FrictionlessComponents`, `Frictional`, `FrictionalPureSlip` (+ optional suffix `WithNormalUpdate`) | Scalar-multiplier ALM (thesis §4.3.3.2.1), vector-multiplier ALM (§4.3.3.2.2, condensable, needed by the `MixedULMLinearSolver`), Coulomb frictional ALM (§4.3.4). `PureSlip` forces every active node into the slip state (`pure_slip` of the criteria). |
| `not_normal_update_frictional` | bool, `false` | Frictional problems update the normals every iteration by default (`WithNormalUpdate` is appended to the contact type); set to `true` to keep the normals of the beginning of the step. |
| `normal_variation` | `no_derivatives_computation` (default), `elemental_derivatives`, `nodal_elemental_derivatives`, `no_derivatives_computation_with_normal_update` (upper-case spellings accepted) | How the slave normals enter the linearisation (`CONSIDER_NORMAL_VARIATION`, enum `NormalDerivativesComputation`). `nodal_elemental_derivatives` selects the `NV` conditions, whose generated tangent includes $$\Delta\mathbf{n}$$ (thesis §4.6.1.4, §4.6.2.4); the others only decide whether the paired normal is refreshed each iteration. |
| `frictional_law` | `Coulomb` (default) | Name of the frictional law. Accepted for future use: the conditions currently implement Coulomb friction directly and Tresca is not wired (see [Frictional laws and MPC constraint](../Implementation/Frictional_Laws_And_MPC_Constraint.html)). |
| `alternative_formulations.axisymmetric` | bool, `false` | Uses the `…Axisym…` conditions (2D only; integrates with $$2\pi r / t$$, `THICKNESS` from the properties). Not available for `FrictionlessComponents`. |
| `interval` | `[0.0, "End"]` | Time interval in which the process is active; outside it the conditions are deactivated. |

### Augmented Lagrangian parameters (`advance_ALM_parameters`)

The ALM functional (thesis eq. 4.9) contains the scale factor $$k$$ and the penalty $$\varepsilon$$; the solution does not depend on them but the conditioning and the convergence do (thesis §4.3.3.3, Tables 4.1–4.2). By default both are computed automatically from the interface stiffness and mesh size (thesis eq. 4.11):

<p align="center">$$ \varepsilon = k \approx \text{stiffness\_factor} \cdot \frac{E_{mean}}{h_{mean}}, \qquad \varepsilon \leftarrow \text{penalty\_scale\_factor} \cdot \varepsilon $$</p>

| Key | Default | Meaning |
|---|---|---|
| `manual_ALM` | `false` | `false`: `ALMVariablesCalculationProcess` computes `INITIAL_PENALTY` and `SCALE_FACTOR` from `YOUNG_MODULUS` and `NODAL_H` of the interface. `true`: the values of `penalty` and `scale_factor` are used verbatim. |
| `stiffness_factor` | `1.0` | Multiplier of $$E_{mean}/h_{mean}$$ for both $$\varepsilon$$ and $$k$$ (the thesis suggests values of order 10). |
| `penalty_scale_factor` | `1.0` | Additional multiplier applied to $$\varepsilon$$ only. |
| `use_scale_factor` | `true` | `false` sets $$k = 1$$ after the automatic computation. |
| `penalty`, `scale_factor` | `1.0e-12`, `1.0` | Manual values (only with `manual_ALM`). A vanishing penalty is replaced by 1. |
| `adapt_penalty` | `false` | Adapted augmented Lagrangian (`AALMAdaptPenaltyValueProcess`, thesis Algorithm 7): the nodal penalty is rescaled from the gap evolution every iteration. |
| `max_gap_factor` | `1.0e-3` | Reference gap of the adaptation as a fraction of `NODAL_H` (`MAX_GAP_FACTOR`). |

The values used are printed at start-up (`SCALE_FACTOR: …`, `INITIAL_PENALTY: …`).

### Frictional parameters

| Key | Default | Meaning |
|---|---|---|
| `tangent_factor` | `2.5e-2` | Ratio between the tangential and the normal penalty: $$\varepsilon_\tau = $$ `tangent_factor` $$\cdot\,\varepsilon$$ (`TANGENT_FACTOR`, thesis §4.3.4.2.3). Penalty processes default to `1.0e-3`, explicit to `1.0e-4`. |
| `slip_threshold` | `2.0e-2` | Hysteresis of the stick/slip decision: a slip node returns to stick only when $$\Vert\bar\lambda_\tau\Vert / (\mu \vert\bar\lambda_n\vert)$$ falls below $$1 - $$ `slip_threshold` (`SLIP_THRESHOLD`, `ActiveSetUtilities::ComputeALMFrictionalActiveSet`). |
| `operator_threshold` | `1.0e-3` | Threshold on $$\Vert\mathbf{D}-\mathbf{D}_{old}\Vert$$ and $$\Vert\mathbf{M}-\mathbf{M}_{old}\Vert$$ that switches the slip increment between the objective and the non-objective measures (`OPERATOR_THRESHOLD`, thesis eqs. 4.65–4.69). |
| `slip_augmentation_coefficient` | `0.0` | Scales the penalty part of the tangential augmented pressure on slip nodes (`SLIP_AUGMENTATION_COEFFICIENT`). |
| `slip_step_reset_frequency` | `1` | Every how many steps the `SLIP` flags are reset (all nodes start the step as stick). `0` never resets; a negative value recomputes the tangents from the `WEIGHTED_SLIP` direction. |
| `friction_coefficients` | `0.0` per interface | $$\mu$$ of each interface. When all coefficients are zero, `AuxiliaryPureSlipCheck` decides whether the problem is run as pure slip. |

### Integration and geometry

| Key | Default | Meaning |
|---|---|---|
| `integration_order` | `2` | Gauss order per integration cell (`INTEGRATION_ORDER_CONTACT`, 1–5); see [Mortar integration](../Theory/Mortar_Integration_And_Dual_Lagrange_Multipliers.html). |
| `consider_tessellation` | `false` | Tessellate warped quadrilaterals before the segmentation (`CONSIDER_TESSELLATION`). Mesh tying defaults to `true`. |
| `zero_tolerance_factor` | `1.0` | Multiplies the machine epsilon used as geometric tolerance in the segmentation (`ZERO_TOLERANCE_FACTOR`). |
| `normal_check_proportion` | `0.1` | Offset (fraction of the element size) used by `NormalCheckProcess` to detect inverted condition normals. |
| `clear_inactive_for_post` | `true` | Zero the augmented pressures of inactive nodes before writing output. |

### `search_parameters`

The search is documented in [Search pipeline and bounding volumes](../Contact_Search/Search_Pipeline_And_Bounding_Volumes.html) and [Gap computation](../Contact_Search/Gap_Computation.html); the table only maps the keys.

| Key | Default | Meaning |
|---|---|---|
| `type_search` | `in_radius_with_obb` | Broad-phase structure: `in_radius`, `in_box`, `in_radius_with_obb`, `in_box_with_obb` (KD-tree, optionally followed by an OBB check), `octree_with_obb`. |
| `simple_search` | `false` | `SimpleContactSearchProcess` (gap-only activation) instead of `AdvancedContactSearchProcess`. |
| `adapt_search` | `false` | Multiplies `search_factor` and `active_check_factor` by the relative mesh-size factor of the interface (`ContactUtilities::CalculateRelativeSizeMesh`). |
| `search_factor` | `3.5` | Search radius / box size as a multiple of the condition size (`NODAL_H`). |
| `active_check_factor` | `0.01` | A slave node is predicted active when its gap is below `active_check_factor` $$\cdot$$ `NODAL_H` (`ACTIVE_CHECK_FACTOR`). |
| `max_number_results`, `bucket_size` | `1000`, `4` | KD-tree parameters (`allocation_size`, `bucket_size`). |
| `dynamic_search` | `false` | Search on the positions predicted with the velocity (needs `VELOCITY`). |
| `static_check_movement` | `false` | Static variant of the movement check. |
| `database_step_update` | `1` | Search frequency in steps (the pairs are kept in between). |
| `normal_orientation_threshold` | `1.0e-1` | Pairs whose unit normals differ by less than this norm are discarded (parallel facets). |
| `consider_gap_threshold` | `false` | Reject pairs whose gap exceeds `DISTANCE_THRESHOLD`. |
| `predict_correct_lagrange_multiplier` | `false` | Predict / correct the multipliers from a gap–pressure regression (advanced search). |
| `check_gap` | `check_mapping` | Gap check mode: `no_check`, `direct_check`, `check_mapping` / `mapping_check` (mortar mapper, `NormalGapProcess`). |
| `debug_mode` | `false` | GiD dumps of the pairs and flags, integration-area report, total load / reaction / contact force print. |
| `octree_search_parameters` | see defaults | Oriented-bounding-box options: `OBB_intersection_type` (`SeparatingAxisTheorem` or `Direct`), `bounding_box_factor`, `build_from_bounding_box`, `lower/higher_bounding_box_coefficient`, `debug_obb`. |

### `advance_explicit_parameters`

Used only by `explicit_penalty_contact_process`: `manual_max_gap_theshold` (sic) with `max_gap_threshold` sets `MAX_GAP_THRESHOLD` (the gap at which the dynamic factor rescales the penalty, `ComputeDynamicFactorProcess`); otherwise the mean `NODAL_H` is used. The remaining keys (`automatic_gap_factor`, `max_gap_factor`, `logistic_exponent_factor`) are accepted by the defaults but not read by the current implementation.

## Condition names created

The process chooses the condition **stem**; the search appends the geometry suffix (`2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N`) of every pair.

| Process | `contact_type` | `normal_variation` = `nodal_elemental_derivatives` | `axisymmetric` | Condition stem |
|---|---|---|---|---|
| ALM | `Frictionless` | no / yes | no | `ALMFrictionlessMortarContact` / `ALMNVFrictionlessMortarContact` |
| ALM | `Frictionless` | no / yes | yes | `ALMFrictionlessAxisymMortarContact` / `ALMNVFrictionlessAxisymMortarContact` |
| ALM | `FrictionlessComponents` | no / yes | – | `ALMFrictionlessComponentsMortarContact` / `ALMNVFrictionlessComponentsMortarContact` |
| ALM | `Frictional*` | no / yes | no | `ALMFrictionalMortarContact` / `ALMNVFrictionalMortarContact` |
| ALM | `Frictional*` | no / yes | yes | `ALMFrictionalAxisymMortarContact` / `ALMNVFrictionalAxisymMortarContact` |
| Penalty (implicit or explicit) | `Frictionless` | no / yes | no / yes | `Penalty[NV]Frictionless[Axisym]MortarContact` |
| Penalty (implicit or explicit) | `Frictional*` | no / yes | no / yes | `Penalty[NV]Frictional[Axisym]MortarContact` |
| MPC | any | – | – | `MPCMortarContact` |
| Mesh tying | – | – | – | `MeshTyingMortar` |

## `penalty_contact_process` and `explicit_penalty_contact_process`

Same keys as the ALM process (the class derives from `ALMContactProcess`) with these differences:

| Key | ALM | Penalty | Explicit penalty |
|---|---|---|---|
| `tangent_factor` | `2.5e-2` | `1.0e-3` | `1.0e-4` |
| `advance_ALM_parameters.penalty` (manual) | `1.0e-12` | `1.0e16` | `1.0e16` |
| `advance_ALM_parameters.max_gap_factor` | `1.0e-3` | `5.0e-4` | `1.0e-3` |
| `search_parameters.type_search` | `in_radius_with_obb` | `in_radius_with_obb` | `octree_with_obb` |
| Multiplier DoFs | yes | none | none |
| `advance_explicit_parameters` | ignored | ignored | sets `MAX_GAP_THRESHOLD`; the process also creates its own `ComputeDynamicFactorProcess` |

The penalty formulation never satisfies the constraint exactly (thesis §4.3.3.2.1.2): with the automatic $$\varepsilon \approx E/h$$ the penetration is of the order of the strain times the element size. Increase `stiffness_factor` for a stiffer contact at the price of conditioning.

## `mpc_contact_process`

Keys shared with the ALM process: `model_part_name`, `contact_model_part`, `assume_master_slave`, `contact_property_ids`, `friction_coefficients`, `contact_type` (`Frictionless` / `Frictional`), `not_normal_update_frictional`, `interval`, `normal_variation`, `frictional_law`, `integration_order`, `consider_tessellation`, `normal_check_proportion`, `clear_inactive_for_post`, `search_parameters`. Specific keys:

| Key | Default | Meaning |
|---|---|---|
| `tangent_factor` | `1.0e-1` | Tangential factor of the frictional constraint update. |
| `zero_tolerance_factor` | `1.0e2` | Tolerance factor (looser than in the mortar conditions). |
| `reaction_check_stiffness_factor` | `1.0e-10` | A tied/contacting node is released when the mapped reaction indicates traction larger than this factor times `YOUNG_MODULUS` (`REACTION_CHECK_STIFFNESS_FACTOR`, `MPCContactCriteria`). |
| `update_condition_relation_step` | `false` | Recompute the constraint relation every step instead of keeping it while the pair survives. |

The solver must carry `mpc_contact_settings` instead of `contact_settings`. Details: [Frictional laws and MPC constraint](../Implementation/Frictional_Laws_And_MPC_Constraint.html).

## `mesh_tying_process`

Documented with its full default block in [Mesh tying](../Theory/Mesh_Tying.html). Specific keys: `variable_name` (tied variable, default `DISPLACEMENT`; scalar → `ScalarMeshTying`, vector → `ComponentsMeshTying`), `consider_static_condensation` (parent elements assigned for condensation), `scale_factor_parameters` (`manual_scale_factor`, `stiffness_factor`, `scale_factor`), `consider_tessellation = true` and `database_step_update = 999999999` by default (the pairing is frozen after the first step).

## Which formulation should I use?

| Situation | Recommendation |
|---|---|
| General frictionless contact, moderate sliding | `alm_contact_process` with `Frictionless` (scalar multiplier, smallest system, block builder). |
| Frictionless contact with many active nodes and iterative solvers | `FrictionlessComponents` + the default `use_mixed_ulm_solver` (multipliers condensed away). |
| Friction | `alm_contact_process` with `Frictional` (or `FrictionalPureSlip` when stick states are irrelevant); expect the buffer size to become 3. |
| Explicit dynamics, impacts | `explicit_penalty_contact_process`. |
| Displacement-only system required, moderate accuracy on non-matching meshes | `mpc_contact_process`. |
| Gluing non-matching meshes | `mesh_tying_process`. |
| Axisymmetric 2D problems | ALM or penalty with `alternative_formulations.axisymmetric = true`. |

See also [Tips, troubleshooting and limitations](Tips_Troubleshooting_And_Limitations.html) for parameter calibration and the [tutorial](Tutorial_Hertz_2D.html) for complete input files.
