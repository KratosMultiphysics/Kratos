# Strategies, convergence criteria and builder-and-solvers

Header-only templates that drive the semi-smooth Newton solution of the contact problem. They are instantiated only in the Python bindings (`custom_python/add_custom_strategies_to_python.cpp`), so this folder contributes no object files to `KratosContactStructuralMechanicsCore`.

![One time step of the contact solution loop](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/images/csma_solution_loop.svg)

## Strategies (`custom_strategies/`)

| Class | Registered name | Adds over the base strategy |
|---|---|---|
| `ResidualBasedNewtonRaphsonContactStrategy` | `newton_raphson_contact_strategy` | Contact-aware `Predict()` (resets `WEIGHTED_GAP`/`WEIGHTED_SLIP`, recomputes the gap through the explicit contribution, moves the nodes); outer loop of the *simplified* semi-smooth Newton when the model part is not `INTERACTION` (`inner_loop_iterations`); adaptive time stepping on failure (`adaptative_strategy`, `split_factor`, `max_number_splits`); inverted-element guard; drives the Python process lists through `ProcessFactoryUtility`. |
| `LineSearchContactStrategy` | `line_search_contact_strategy` | Split line search: the increment is separated into displacement and Lagrange-multiplier parts and a parabola is fitted to each block (`ComputeSplitDx`, `ComputeMixedResidual`, `ComputeParabola`). |
| `ResidualBasedNewtonRaphsonMPCContactStrategy` | `newton_raphson_mpc_contact_strategy` | Route based on multipoint constraints: internal `MPCContactCriteria`, `ComputeNodalWeights()` every iteration, one `BuildAndSolve` in `Predict()` to check the active set, `update_each_nl_iteration` (DoF set rebuilt every iteration), `enforce_ntn` accepted but not implemented (`EnforcingNTN()` is commented out). |

Defaults (merged with the base Newton–Raphson defaults):

```json
{
    "name"                  : "newton_raphson_contact_strategy",
    "adaptative_strategy"   : false,
    "split_factor"          : 10.0,
    "max_number_splits"     : 3,
    "inner_loop_iterations" : 5
}
```

```json
{
    "name"                     : "newton_raphson_mpc_contact_strategy",
    "inner_loop_iterations"    : 5,
    "update_each_nl_iteration" : false,
    "enforce_ntn"              : false
}
```

The two semi-smooth Newton modes are selected in the solver settings: `simplified_semi_smooth_newton = false` (default) sets the `INTERACTION` flag on the computing model part and solves displacements, multipliers and active set in one Newton loop; `true` freezes the sets during an inner Newton loop and re-checks them up to `inner_loop_iterations` times.

## Convergence criteria (`custom_convergencecriterias/`)

All derive from `ConvergenceCriteria`, expose `GetDefaultParameters()`, `Name()` and `Create(Parameters)`, and are combined by `python_scripts/contact_convergence_criteria_factory.py` into a `MortarAndConvergenceCriteria` = *user criterion* **and** *mortar (active set) criterion*.

| Header | Class (registered name) | Checks |
|---|---|---|
| `base_mortar_criteria.h` | `BaseMortarConvergenceCriteria` (`base_mortar_criteria`) | Base of the mortar criteria: `PreCriteria` updates normals/tangents, recomputes the weighted gap, runs `ComputeDynamicFactorProcess` and `AALMAdaptPenaltyValueProcess`; `PostCriteria` stores the gap history and writes the GiD debug output (`gidio_debug`). |
| `alm_frictionless_mortar_criteria.h` | `ALMFrictionlessMortarConvergenceCriteria` | Active set of the scalar ALM (`ActiveSetUtilities::ComputeALMFrictionlessActiveSet`). |
| `alm_frictionless_components_mortar_criteria.h` | `ALMFrictionlessComponentsMortarConvergenceCriteria` | Active set of the vector (components) ALM. |
| `alm_frictional_mortar_criteria.h` | `ALMFrictionalMortarConvergenceCriteria` | Active set and stick/slip set of the frictional ALM (`pure_slip` option). |
| `penalty_frictionless_mortar_criteria.h`, `penalty_frictional_mortar_criteria.h` | `PenaltyFrictionless…`, `PenaltyFrictional…MortarConvergenceCriteria` | Same for the penalty formulations. |
| `mesh_tying_mortar_criteria.h` | `MeshTyingMortarConvergenceCriteria` | No active set; table hook only (always converged). |
| `mpc_contact_criteria.h` | `MPCContactCriteria` | Active set of the MPC route from the reactions mapped master → slave with the mortar mapper (`REACTION_CHECK_STIFFNESS_FACTOR`). |
| `mortar_and_criteria.h` | `MortarAndConvergenceCriteria` (`mortar_and_criteria`) | `AndCriteria` specialisation that prints the convergence table and, optionally, the condition number of the system. |
| `displacement_contact_criteria.h`, `displacement_residual_contact_criteria.h` | `DisplacementContactCriteria`, `DisplacementResidualContactCriteria` | Displacement (+ rotation) increment or residual — used by the penalty formulations (no multipliers). |
| `displacement_lagrangemultiplier_{,residual_,mixed_}contact_criteria.h` | `DisplacementLagrangeMultiplier{,Residual,Mixed}ContactCriteria` | Displacement and Lagrange multiplier increment / residual / mixed (displacement residual + LM increment) with separate tolerances (`contact_*_tolerance`). |
| `displacement_lagrangemultiplier_{,residual_,mixed_}frictional_contact_criteria.h` | `…FrictionalContactCriteria` | Same, with the contact block split into stick and slip nodes (`frictional_stick_*`, `frictional_slip_*`, `ratio_normal_tangent_threshold`). The residual variant is registered as `displacement_lagrangemultiplier_ressidual_frictional_contact_criteria` (typo kept in the source). |
| `contact_error_mesh_criteria.h` | `ContactErrorMeshCriteria` | Discretisation error (`ContactSPRErrorProcess`) against `error_mesh_tolerance`; "converged" means no remeshing needed (adaptive remeshing, needs the `MeshingApplication`). |

Default tolerances are `1.0e-4` (relative) and `1.0e-9` (absolute) for every block; `ensure_contact` fails the step if no node is active; `print_convergence_criterion` prints the table of every criterion.

## Builder-and-solvers (`custom_builder_and_solvers/`)

| Class | Registered name | Purpose |
|---|---|---|
| `ContactResidualBasedBlockBuilderAndSolver` | `contact_block_builder_and_solver` | Block B&S that fixes the Lagrange-multiplier DoFs of `ISOLATED` nodes (nodes whose pairs have all been deactivated) before applying the Dirichlet conditions and frees them afterwards, so the multiplier block never becomes singular. |
| `ContactResidualBasedEliminationBuilderAndSolver` | `contact_residual_elimination_builder_and_solver` | Elimination B&S that also fixes the multiplier DoF of a slave node whose displacement is fixed. |
| `ContactResidualBasedEliminationBuilderAndSolverWithConstraints` | `contact_residual_elimination_builder_and_solver_with_constraints` | Same with master–slave constraints (the MPC route): `SetUpSystemWithConstraints`, DoF set including the constraint DoFs. |

The condensed solution of the dual-multiplier system is not done here but in [`../custom_linear_solvers/`](../custom_linear_solvers/README.md) (`MixedULMLinearSolver`).

## Full documentation

- [Strategies and convergence criteria](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Strategies_And_Convergence_Criteria.html) · [source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Strategies_And_Convergence_Criteria.md)
- [Builder-and-solvers and linear solvers](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Builder_And_Solvers_And_Linear_Solvers.html) · [source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Builder_And_Solvers_And_Linear_Solvers.md)
- [Solver settings reference](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Solver_Settings_Reference.html) · [source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Usage/Solver_Settings_Reference.md)
- Theory of the semi-smooth Newton method: [frictionless](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Frictionless_Contact.md) (thesis Algorithm 2) and [frictional](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Frictional_Contact.md) (Algorithm 3) contact
