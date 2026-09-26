# Python scripts

The Python layer of the application: the contact solvers (thin extensions of the `StructuralMechanicsApplication` solvers), the processes that set up the contact interface and drive the search, and the helpers that build criteria, strategies and linear solvers from the JSON settings.

## How a contact simulation is wired

`StructuralMechanicsAnalysis` is reused unchanged. The generic wrapper `StructuralMechanicsApplication/python_scripts/python_solvers_wrapper_structural.py` picks a solver of this folder as soon as `solver_settings` contains a `contact_settings` block (→ `contact_<static|implicit_dynamic|explicit_dynamic>_solver`) or an `mpc_contact_settings` block (→ `mpc_contact_<static|implicit_dynamic>_solver`). The contact itself is set up by one process of `processes.contact_process_list`, which also writes `contact_settings.mortar_type` into the solver settings.

```json
"solver_settings" : {
    "solver_type"           : "Static",
    "contact_settings"      : { "mortar_type" : "ALMContactFrictionless" },
    "convergence_criterion" : "contact_residual_criterion"
},
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

## Solvers

| File | Class | Notes |
|---|---|---|
| `contact_structural_mechanics_static_solver.py` | `ContactStaticMechanicalSolver` | Static contact (Newton–Raphson / line search / arc length contact strategies). |
| `contact_structural_mechanics_implicit_dynamic_solver.py` | `ContactImplicitMechanicalSolver` | Implicit dynamics (Newmark / Bossak), `compute_dynamic_factor`. |
| `contact_structural_mechanics_explicit_dynamic_solver.py` | `ContactExplicitMechanicalSolver` | Explicit (central difference) with penalty contact; `delta_time_factor_for_contact`. |
| `mpc_contact_structural_mechanics_static_solver.py`, `mpc_contact_structural_mechanics_implicit_dynamic_solver.py` | `MPCContactStaticSolver`, `MPCContactImplicitMechanicalSolver` | Multipoint-constraint route (`mpc_contact_settings`). |
| `adaptive_remeshing/*.py`, `python_solvers_wrapper_adaptative_remeshing_contact_structural.py` | `AdaptativeRemeshingContact*` | Solvers, analysis stage and utilities for error-driven remeshing with MMG (needs the `MeshingApplication`). |
| `auxiliary_methods_solvers.py` | functions | Shared implementation: default `contact_settings` (`AuxiliaryContactSettings`), variables and DoFs per `mortar_type` (`AuxiliaryAddVariables`, `AuxiliaryAddDofs`), strategy and linear-solver creation (`AuxiliaryNewton`, `AuxiliaryLineSearch`, `AuxiliaryCreateLinearSolver` with the `MixedULMLinearSolver` wrapping), buffer-size rule (frictional ⇒ ≥ 3). |
| `contact_convergence_criteria_factory.py` | `ContactConvergenceCriteriaFactory` | Builds *user criterion* AND *mortar criterion* (`MortarAndConvergenceCriteria`) from `convergence_criterion` and `mortar_type`. |

`mortar_type` values: `ALMContactFrictionless`, `ALMContactFrictionlessComponents`, `ALMContactFrictional` (+ `PureSlip`), `PenaltyContactFrictionless`, `PenaltyContactFrictional` (+ `PureSlip`), `ScalarMeshTying`, `ComponentsMeshTying`. They decide the nodal variables and DoFs, the active-set criterion and whether the `MixedULMLinearSolver` is inserted (vector multipliers only).

## Processes

| File | Class | Condition stem | Purpose |
|---|---|---|---|
| `search_base_process.py` | `SearchBaseProcess` | – | Base: builds the `Contact` / `ContactSubN` / `MasterSubModelPartN` / `SlaveSubModelPartN` sub-model-parts, runs `InterfacePreprocessCondition`, `NormalCheckProcess`, `MasterSlaveProcess`, creates one `ContactSearchProcess` per pair and re-runs it every `database_step_update` steps. |
| `alm_contact_process.py` | `ALMContactProcess` | `ALM[NV]Frictionless[Components][Axisym]MortarContact`, `ALM[NV]Frictional[Axisym]MortarContact` | Augmented Lagrangian contact (`contact_type` Frictionless / FrictionlessComponents / Frictional); automatic $\varepsilon$, $k$ (`advance_ALM_parameters`), frictional parameters, normal variation, axisymmetry. |
| `penalty_contact_process.py` | `PenaltyContactProcess` | `Penalty[NV]Frictionless|Frictional[Axisym]MortarContact` | Penalty contact (same keys, penalty defaults). |
| `explicit_penalty_contact_process.py` | `ExplicitPenaltyContactProcess` | idem | Penalty contact for the explicit solver (`octree_with_obb` search by default). |
| `mpc_contact_process.py` | `MPCContactProcess` | `MPCMortarContact` | Multipoint-constraint contact (`reaction_check_stiffness_factor`, `update_condition_relation_step`). |
| `mesh_tying_process.py` | `MeshTyingProcess` | `MeshTyingMortar` | Mortar mesh tying of a `variable_name` (scalar or vector), optional static condensation. |
| `contact_remesh_mmg_process.py` | `ContactRemeshMmgProcess` | – | MMG remeshing adapted to contact (metric on stress / contact pressure / strain energy; needs `MeshingApplication`). |
| `basic_mapping_process.py` | `BasicMappingProcess` | – | Thin wrapper over `SimpleMortarMapperProcess` with interval support. |
| `replace_properties_process.py` | `ReplacePropertiesProcess` | – | Reloads materials from a JSON at a given interval. Note: its default-settings string contains `"reinitialize_entities" : false.` (period instead of comma), which fails to parse when the defaults are used. |
| `custom_sympy_fe_utilities.py` | functions | – | Symbolic helpers for the code generators of `../automatic_differentiation/` (sympy 1.2). |

## Full documentation

- [Getting started](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/General/Getting_Started.html) · [source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/General/Getting_Started.md)
- [Solver settings reference](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Solver_Settings_Reference.md) · [Contact process settings reference](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Contact_Process_Settings_Reference.md) · [Tutorial](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Tutorial_Hertz_2D.md) · [Output and post-processing](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Output_And_Postprocessing.md) · [Tips, troubleshooting and limitations](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Tips_Troubleshooting_And_Limitations.md)
- [Adaptive remeshing](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Examples/Adaptive_Remeshing.md)
