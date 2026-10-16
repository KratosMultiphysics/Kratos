---
title: Strategies and Convergence Criteria
keywords: contact, strategy, Newton-Raphson, semi-smooth Newton, active set, convergence criteria, line search, MPC, INTERACTION, adaptative time step, Kratos
tags: [contact, implementation, strategies, convergence, ContactStructuralMechanicsApplication]
sidebar: contact_structural_mechanics_application
summary: The three solving strategies of the ContactStructuralMechanicsApplication (contact Newton-Raphson, contact line search and MPC contact), the 18 convergence criteria that drive the active set and measure convergence, how the Python factory composes them from the solver settings, and how all of it maps to the semi-smooth Newton algorithms of the thesis.
---

> **Sources.** Thesis §4.3.3.5–4.3.3.6 (pp. 111–112, Algorithm 2, eqs. 4.41–4.44), §4.3.4.4–4.3.4.5 (pp. 122–123, Algorithm 3, eqs. 4.73–4.79), Algorithm 1 (line search); code: `custom_strategies/custom_strategies/residualbased_newton_raphson_contact_strategy.h`, `line_search_contact_strategy.h`, `residualbased_newton_raphson_mpc_contact_strategy.h`, `custom_strategies/custom_convergencecriterias/*.h` (18 headers), `custom_utilities/active_set_utilities.cpp`, `custom_python/add_custom_strategies_to_python.cpp`, `python_scripts/contact_convergence_criteria_factory.py`, `python_scripts/auxiliary_methods_solvers.py`, `python_scripts/contact_structural_mechanics_static_solver.py`, `python_scripts/mpc_contact_structural_mechanics_static_solver.py`.

The contact problem is solved with a *semi-smooth Newton* method: the non-linearity of the contact constraints (which nodes are in contact, which of them stick and which slip) is handled by a primal-dual active-set update that runs together with the standard Newton–Raphson iterations of the structural problem. In the ContactStructuralMechanicsApplication this is split into two families of header-only classes living in `custom_strategies/`:

- the **strategies** (`custom_strategies/custom_strategies/`) decide *how the Newton loop is organized*: whether the active set is updated inside every iteration or in an outer loop, whether the time step is split when the iterations fail, how the solution is predicted, and how the line search treats displacements and Lagrange multipliers;
- the **convergence criteria** (`custom_strategies/custom_convergencecriterias/`) decide *when the loop stops*: they measure the displacement / residual / multiplier norms, but also update the mortar quantities (normals, tangents, weighted gap) and perform the active-set check itself.

All of them are templates on `TSparseSpace, TDenseSpace[, TLinearSolver]` instantiated only for the Ublas spaces in [`add_custom_strategies_to_python.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_python/add_custom_strategies_to_python.cpp); no time-integration *scheme* is defined by the application (the `SCHEME CLASSES` block of that file is empty), the schemes of the StructuralMechanicsApplication are reused. How the pieces are wired from `ProjectParameters.json` is summarized at the end of this page and detailed in the [Solver settings reference](../Usage/Solver_Settings_Reference.html); the mathematics is in [Frictionless contact](../Theory/Frictionless_Contact.html) and [Frictional contact](../Theory/Frictional_Contact.html).

<p align="center"><img src="images/csma_solution_loop.svg" alt="Contact solution loop: search, predict, Newton iterations with PreCriteria / PostCriteria and active-set update, finalization" width="1000"/></p>
<p align="center"><em>Figure: one time step of the contact solution loop and where the strategy and the convergence criteria intervene.</em></p>

## Strategies

| Class | File | `Name()` | Python name | Base class | Created by |
|---|---|---|---|---|---|
| `ResidualBasedNewtonRaphsonContactStrategy` | `residualbased_newton_raphson_contact_strategy.h` | `newton_raphson_contact_strategy` | `ResidualBasedNewtonRaphsonContactStrategy` | `ResidualBasedNewtonRaphsonStrategy` | `auxiliary_methods_solvers.AuxiliaryNewton` when `solving_strategy_settings.type` is `newton_raphson` and `contact_settings.mortar_type` is not empty |
| `LineSearchContactStrategy` | `line_search_contact_strategy.h` | `line_search_contact_strategy` | `LineSearchContactStrategy` | `LineSearchStrategy` | `auxiliary_methods_solvers.AuxiliaryLineSearch` when `solving_strategy_settings.type` is `line_search` |
| `ResidualBasedNewtonRaphsonMPCContactStrategy` | `residualbased_newton_raphson_mpc_contact_strategy.h` | `newton_raphson_mpc_contact_strategy` | `ResidualBasedNewtonRaphsonMPCContactStrategy` | `ResidualBasedNewtonRaphsonStrategy` | `auxiliary_methods_solvers.AuxiliaryMPCNewton` from the `MPCContactStaticSolver` / `MPCContactImplicitMechanicalSolver` |

The `Name()` string is the value of the `"name"` key of `GetDefaultParameters()`; it is the identifier the Kratos strategy factories use through `Create(rModelPart, ThisParameters)`. The `arc_length` option of `solving_strategy_settings.type` falls back to the structural `_create_arc_length_strategy`, which is not contact-aware.

### `ResidualBasedNewtonRaphsonContactStrategy`

[Source](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_strategies/custom_strategies/residualbased_newton_raphson_contact_strategy.h). It derives from the core `ResidualBasedNewtonRaphsonStrategy` and keeps all of its constructors (with or without a builder-and-solver, with `Parameters`), adding two optional trailing arguments `pMyProcesses` and `pPostProcesses` of type `ProcessFactoryUtility::Pointer` (`ProcessesListType`). Its default parameters are merged recursively with the base-class defaults:

```json
{
    "name"                                : "newton_raphson_contact_strategy",
    "adaptative_strategy"                 : false,
    "split_factor"                        : 10.0,
    "max_number_splits"                   : 3,
    "inner_loop_iterations"               : 5
}
```

The Python solvers copy exactly these four keys (`adaptative_strategy`, `split_factor`, `max_number_splits`, `inner_loop_iterations`) from `contact_settings` into the strategy parameters (`AuxiliaryNewton`, `auxiliary_methods_solvers.py`). What the class changes with respect to the base strategy:

**Contact-aware `Predict()`** (`residualbased_newton_raphson_contact_strategy.h:313-348`). The base predictor is *not* called (the call is commented out with the remark "May cause problems in dynamics"). Instead, if the `Contact` sub-model-part nodes carry `WEIGHTED_GAP`, the strategy

1. sets `WEIGHTED_GAP` to zero on the `Contact` nodes, and `WEIGHTED_SLIP` to zero too when the model part `Is(SLIP)` (frictional problem);
2. calls `ContactUtilities::ComputeExplicitContributionConditions` on `ComputingContact`, so that every pair condition integrates its mortar operators and accumulates the current weighted gap $$\tilde{g}_n$$ on the slave nodes;
3. moves the nodal coordinates by the displacement increment: `DISPLACEMENT` at `STEP == 1`, `DISPLACEMENT - DISPLACEMENT(1)` afterwards.

No prediction of the Lagrange multipliers is made here (the corresponding block is commented out); the `predict_correct_lagrange_multiplier` option of the search process takes care of it.

**Two ways of running the semi-smooth Newton loop** (`SolveSolutionStep`, lines 465-515). The choice is made by the `INTERACTION` flag of the computing model part, which `ContactStaticMechanicalSolver.Initialize` sets to `true` unless `contact_settings.simplified_semi_smooth_newton` is `true`:

| Computing model part | Loop | Active-set check |
|---|---|---|
| `Is(INTERACTION)` (default) | One call to `BaseSolveSolutionStep()`; `INNER_LOOP_ITERATION = 1` | Inside every Newton iteration, by the mortar criterion's `PostCriteria` (full semi-smooth Newton, thesis §4.3.3.6) |
| `IsNot(INTERACTION)` (`simplified_semi_smooth_newton`) | Up to `inner_loop_iterations` outer passes; each pass resets `NL_ITERATION_NUMBER = 1`, stores the pass counter in `INNER_LOOP_ITERATION`, runs a complete `BaseSolveSolutionStep()` with a *frozen* active set, then calls `mpConvergenceCriteria->PostCriteria` (echo temporarily silenced) to update the set and decide whether another pass is needed | Only at the end of each inner pass; `ActiveSetUtilities` only re-evaluates the set when `rModelPart.Is(INTERACTION) || NL_ITERATION_NUMBER == 1` |

**Own Newton loop `BaseSolveSolutionStep()`** (lines 615-790). It is a copy of the core Newton–Raphson loop instrumented for contact: `PreCriteria` is evaluated before the first build, `InitializeNonLinearIteration` / `FinalizeNonLinearIteration` of the criteria are called around every iteration, the RHS is rebuilt before `PostCriteria` when `GetActualizeRHSflag()` is set (the factory sets `SetActualizeRHSFlag(True)`), and the geometry is checked for inversion when `adaptative_strategy` is on. The echo of the convergence criterion is silenced during `InitializeSolutionStep`, and the `mFinalizeWasPerformed` guard ensures that `FinalizeSolutionStep` runs once even when the step has been split.

**Adaptive time-step splitting** (`AdaptativeStep()`, lines 796-925; `SplitTimeStep`, `UnMoveMesh`, `CoutSplittingTime`). Only when `adaptative_strategy` is `true` and the step did not converge. The strategy restores `TIME` to the beginning of the step, divides `DELTA_TIME` by `split_factor`, and marches through the original interval in sub-steps; the process is repeated (dividing again) up to `max_number_splits` times. Each sub-step increments `STEP`, overwrites/clones the solution-step data in the nodal buffer and the `ProcessInfo`, and runs the full sequence `InitializeSolutionStep → Predict → SolveSolutionStep → FinalizeSolutionStep`. Around each sub-step the two `ProcessFactoryUtility` lists are driven explicitly: `ExecuteInitializeSolutionStep`, `ExecuteFinalizeSolutionStep`, `ExecuteBeforeOutputStep`, `PrintOutput()` (post processes) and `ExecuteAfterOutputStep`. This is why `ContactStaticMechanicalSolver.AddProcessesList(processes_list)` and `AddPostProcess(post_process)` exist: they wrap the Python process lists in a `ProcessFactoryUtility` so that boundary conditions depending on time and the output are still applied for the sub-steps. If the lists are missing, the strategy warns that the adaptive strategy "will be USELESS" (echo level greater than 0). After the last split the original `DELTA_TIME` is restored; if the splits are exhausted, `MaxIterationsAndSplitsExceeded()` prints a warning box.

**Inverted-element guard** (`CheckGeometryInverted()`, lines 928-964). Also only with `adaptative_strategy`. Before the first solve of the step and after every `UpdateDatabase`, the strategy loops over the elements (serially) and returns `true` if any `DeterminantOfJacobian(0)` is negative or if any `DEFORMATION_GRADIENT` computed on the integration points has a negative determinant. In that case `STEP` is decremented, a `KRATOS_WARNING("Element inverted")` is issued and `BaseSolveSolutionStep` returns `false`, so the adaptive splitting kicks in.

**Extra Python methods.** Besides the usual `Initialize`, `Solve`, `Predict`, `SolveSolutionStep`, … inherited from `SolvingStrategy`, the binding exposes `SetMaxIterationNumber`, `GetMaxIterationNumber`, `SetKeepSystemConstantDuringIterations` and `GetKeepSystemConstantDuringIterations` (the latter controls whether `BuildAndSolve` or `BuildRHSAndSolve` is called in the iterations).

### `LineSearchContactStrategy`

[Source](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_strategies/custom_strategies/line_search_contact_strategy.h). It derives from the core `LineSearchStrategy` (thesis Algorithm 1) and has only one extra default parameter:

```json
{
    "name" : "line_search_contact_strategy"
}
```

The core line search fits one parabola to the norm of the *whole* residual at the relaxation factors 0, ½ and 1. In a mixed displacement / Lagrange-multiplier system the two blocks have very different magnitudes, so the contact version overrides `UpdateDatabase` (lines 350-403) with a **split line search**:

1. `ComputeSplitDx(Dx, DxDisp, DxLM)` separates the solution increment by DoF variable: `DISPLACEMENT_X/Y/Z` go to `DxDisp`, everything else (the multipliers) to `DxLM`.
2. `ComputeMixedResidual(b, normDisp, normLM)` returns *two* residual norms with the same split. The residual is evaluated without update ($$r_o$$), after half of the increment ($$r_h$$) and after the full increment ($$r_f$$), each time rebuilding the RHS with `BuildRHS`.
3. `ComputeParabola` is called once per block with the triple $$(r_f, r_o, r_h)$$ and fits $$y = a x^2 + b x + c$$ with $$c = r_o$$, $$b = 4 r_h - r_f - 3 r_o$$, $$a = 2 r_f - 4 r_h + 2 r_o$$. If $$a \gt 0$$ the minimum $$x = -b/(2a)$$ is taken and clamped to $$[-1, 1]$$; otherwise the full step is kept when $$r_f \lt r_o$$ and the lower bound `Xmin` ($$10^{-3}$$, "should be zero, but otherwise it will stagnate") is used when it is not.
4. A final `UpdateDatabase` applies $$-(1 - X_{disp})\,\Delta\mathbf{u}$$ and $$-(1 - X_{LM})\,\Delta\boldsymbol{\lambda}$$ on top of the already applied full step, so that the two blocks end up with independent relaxation factors.

Note: the signature is `ComputeParabola(double& Xmax, double& Xmin, rf, ro, rh)` but the two calls pass `(XminDisp, XmaxDisp, ...)` and `(XminLM, XmaxLM, ...)`, so the optimum is written into the `Xmin*` variables while the final correction uses `XmaxDisp` / `XmaxLM`, which stay at 1.0 (except in the degenerate branch, where `Xmin` is assigned the value 1.0 of the swapped argument). In practice the final correction is therefore zero and the full Newton step is kept; the split parabola machinery is evaluated but does not alter the update. The member `mRecalculateFactor` is declared but not used by the override. The strategy also exposes `SetMaxIterationNumber`, `GetMaxIterationNumber`, `SetKeepSystemConstantDuringIterations` and `GetKeepSystemConstantDuringIterations` to Python; the Python wrapper `AuxiliaryLineSearch` passes an empty `Parameters` object and ignores the processes lists.

### `ResidualBasedNewtonRaphsonMPCContactStrategy`

[Source](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_strategies/custom_strategies/residualbased_newton_raphson_mpc_contact_strategy.h). This is the strategy of the multi-point-constraint (MPC) contact route, in which the contact constraints are imposed as `ContactMasterSlaveConstraint` objects instead of Lagrange multipliers or penalties (see [Frictional laws and MPC constraint](Frictional_Laws_And_MPC_Constraint.html)). Defaults:

```json
{
    "name"                                : "newton_raphson_mpc_contact_strategy",
    "inner_loop_iterations"               : 5,
    "update_each_nl_iteration"            : false,
    "enforce_ntn"                         : false
}
```

`AuxiliaryMPCNewton` copies `inner_loop_iterations`, `update_each_nl_iteration` and `enforce_ntn` from `mpc_contact_settings` (whose own default for `inner_loop_iterations` is 10). Behavior:

- **Internal active-set criterion.** Every constructor creates its own `MPCContactCriteria` (`mpMPCContactCriteria`), used *independently of the user criterion* for the active-set check. The Python solver additionally wraps the structural criterion and a second `MPCContactCriteria` into a `KM.AndCriteria`, so the mapped-reaction check is performed twice per iteration.
- **`Predict()`** (lines 311-345): calls the base predictor and then performs one full `BuildAndSolve` followed by `mpMPCContactCriteria->PostCriteria(...)`, so that the active set of the step starts from the reactions of a first linear solve rather than from the search-based guess.
- **`ComputeNodalWeights()`** (lines 823-860): on `Initialize`, `InitializeSolutionStep` and at every non-linear iteration, it resets and recomputes two non-historical nodal values on the `Contact` sub-model-part from the `SLAVE` conditions: `NODAL_PAUX` (number of slave conditions sharing the node) and `NODAL_MAUX` (lumped area, `LumpingFactors × DomainSize`). `MPCMortarContactCondition` divides the constraint rows by `NODAL_PAUX` and `MPCContactCriteria` divides the reactions by `NODAL_MAUX` to obtain pressures.
- **Loop structure** (`SolveSolutionStep`, lines 412-457): here the `INTERACTION` flag is read from the *process info* (set by `MPCContactStaticSolver.Initialize` when `mpc_contact_settings.simplified_semi_smooth_newton` is `true`). If set, up to `inner_loop_iterations` passes of `AuxiliarySolveSolutionStep()` are run, each followed by `mpMPCContactCriteria->PostCriteria`; otherwise a single pass. `AuxiliarySolveSolutionStep` is the instrumented Newton loop.
- **`update_each_nl_iteration`**: sets the `INTERACTION` flag on the conditions of `ComputingContact` (`VariableUtils().SetFlag(INTERACTION, update_each_nl_iteration, ...)`), which is what makes `MPCMortarContactCondition::InitializeNonLinearIteration` rebuild its constraint, and inside the Newton loop calls `SetUpDofSet`, `SetUpSystem` and `ResizeAndInitializeVectors` at every iteration, because the set of active constraints (hence the DoF numbering of the builder-and-solver with constraints) changes with the active set.
- **`enforce_ntn`** (node-to-node enforcement): the parameter is accepted and validated, but every call to `EnforcingNTN()` is commented out, the method body itself is commented out, and `ComputeNodalWeights` hard-codes `const bool enforce_ntn = false;`. Note: setting `enforce_ntn` to `true` therefore has no effect in the current code.

## Convergence criteria

Eighteen criteria live in `custom_strategies/custom_convergencecriterias/`. All of them expose `GetDefaultParameters()`, `Name()` and `Create(Parameters)`; the Python names coincide with the class names except for `BaseMortarConvergenceCriteria`, which is not exposed. Every default block below is merged with the base `ConvergenceCriteria` defaults (`RecursivelyAddMissingParameters`).

| Class | Header | `Name()` | Base | What it checks |
|---|---|---|---|---|
| `BaseMortarConvergenceCriteria` | `base_mortar_criteria.h` | `base_mortar_criteria` | `ConvergenceCriteria` | Nothing by itself (always `true`); updates normals, tangents, weighted gap, dynamic factor and adapted penalty |
| `ALMFrictionlessMortarConvergenceCriteria` | `alm_frictionless_mortar_criteria.h` | `alm_frictionless_mortar_criteria` | `BaseMortarConvergenceCriteria` | Active set, scalar ALM |
| `ALMFrictionlessComponentsMortarConvergenceCriteria` | `alm_frictionless_components_mortar_criteria.h` | `alm_frictionless_components_mortar_criteria` | `BaseMortarConvergenceCriteria` | Active set, vector (components) ALM |
| `ALMFrictionalMortarConvergenceCriteria` | `alm_frictional_mortar_criteria.h` | `alm_frictional_mortar_criteria` | `BaseMortarConvergenceCriteria` | Active set and stick/slip set, ALM frictional |
| `PenaltyFrictionlessMortarConvergenceCriteria` | `penalty_frictionless_mortar_criteria.h` | `penalty_frictionless_mortar_criteria` | `BaseMortarConvergenceCriteria` | Active set, penalty frictionless |
| `PenaltyFrictionalMortarConvergenceCriteria` | `penalty_frictional_mortar_criteria.h` | `penalty_frictional_mortar_criteria` | `BaseMortarConvergenceCriteria` | Active set and stick/slip set, penalty frictional |
| `MeshTyingMortarConvergenceCriteria` | `mesh_tying_mortar_criteria.h` | `mesh_tying_mortar_criteria` | `ConvergenceCriteria` | Nothing (table hook, always `true`) |
| `MPCContactCriteria` | `mpc_contact_criteria.h` | `mpc_contact_criteria` | `ConvergenceCriteria` | Active set (and stick/slip) from mapped reactions, MPC route |
| `MortarAndConvergenceCriteria` | `mortar_and_criteria.h` | `mortar_and_criteria` | `And_Criteria` | Both children; prints the iteration table and optionally the condition number |
| `DisplacementContactCriteria` | `displacement_contact_criteria.h` | `displacement_contact_criteria` | `ConvergenceCriteria` | Displacement (and rotation) increment |
| `DisplacementResidualContactCriteria` | `displacement_residual_contact_criteria.h` | `displacement_residual_contact_criteria` | `ConvergenceCriteria` | Displacement (and rotation) residual |
| `DisplacementLagrangeMultiplierContactCriteria` | `displacement_lagrangemultiplier_contact_criteria.h` | `displacement_lagrangemultiplier_contact_criteria` | `ConvergenceCriteria` | Displacement increment and Lagrange-multiplier increment |
| `DisplacementLagrangeMultiplierResidualContactCriteria` | `displacement_lagrangemultiplier_residual_contact_criteria.h` | `displacement_lagrangemultiplier_residual_contact_criteria` | `ConvergenceCriteria` | Displacement residual and Lagrange-multiplier residual |
| `DisplacementLagrangeMultiplierMixedContactCriteria` | `displacement_lagrangemultiplier_mixed_contact_criteria.h` | `displacement_lagrange_multiplier_mixed_contact_criteria` | `ConvergenceCriteria` | Displacement *residual* and Lagrange-multiplier *increment* |
| `DisplacementLagrangeMultiplierFrictionalContactCriteria` | `displacement_lagrangemultiplier_frictional_contact_criteria.h` | `displacement_lagrangemultiplier_frictional_contact_criteria` | `ConvergenceCriteria` | Displacement increment; multiplier increment split into normal, stick-tangent and slip-tangent |
| `DisplacementLagrangeMultiplierResidualFrictionalContactCriteria` | `displacement_lagrangemultiplier_residual_frictional_contact_criteria.h` | `displacement_lagrangemultiplier_ressidual_frictional_contact_criteria` (sic) | `ConvergenceCriteria` | Residual version of the previous one |
| `DisplacementLagrangeMultiplierMixedFrictionalContactCriteria` | `displacement_lagrangemultiplier_mixed_frictional_contact_criteria.h` | `displacement_lagrangemultiplier_mixed_frictional_contact_criteria` | `ConvergenceCriteria` | Mixed version of the previous one |
| `ContactErrorMeshCriteria` | `contact_error_mesh_criteria.h` | `contact_error_mesh_criteria` | `ConvergenceCriteria` | Mesh discretization error (SPR); "converged" means *no remeshing needed* |

Note: the `Name()` of `DisplacementLagrangeMultiplierResidualFrictionalContactCriteria` is misspelled `..._ressidual_...` in the source (`displacement_lagrangemultiplier_residual_frictional_contact_criteria.h:625,653`); any factory-based creation must use the misspelled string.

The following conventions are shared by the whole family:

- **Relative and absolute checks.** Every tolerance comes in a pair `*_relative_tolerance` / `*_absolute_tolerance`; a block converges when *either* the ratio $$\sqrt{\Vert\Delta x\Vert^2/\Vert x\Vert^2}$$ (or the residual ratio with respect to the residual of the first iteration) is below the relative tolerance *or* the absolute value $$\sqrt{\Vert\Delta x\Vert^2}/n_{dof}$$ is below the absolute tolerance. All defaults are $$10^{-4}$$ (relative) and $$10^{-9}$$ (absolute).
- **DoF classification by variable.** The DoF set is partitioned with the variable of each DoF: `DISPLACEMENT_X/Y/Z` form the displacement block, `ROTATION_X/Y/Z` the rotation block (only when the model part has rotation DoFs, local flag `ROTATION_DOF_IS_CONSIDERED`), `LAGRANGE_MULTIPLIER_CONTACT_PRESSURE` or `VECTOR_LAGRANGE_MULTIPLIER_X/Y/Z` the contact block. In the frictional variants the vector multiplier of each slave node is projected on the nodal `NORMAL` to obtain the normal component, and the tangential remainder is assigned to the *stick* or *slip* block according to the `SLIP` flag of the node (or to slip for everybody when `pure_slip` is `true`).
- **`ensure_contact`** (local flag `ENSURE_CONTACT`). When `false` (default) and the norm of the multipliers is zero (no active node), the multiplier block is considered converged; when `true` a zero multiplier norm raises `KRATOS_ERROR` with the message "CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?".
- **`ratio_normal_tangent_threshold`** (frictional variants). A tangential block that fails its own tolerances is still accepted if the ratio between its absolute norm and the absolute norm of the normal block is below this threshold (default $$10^{-4}$$), because the tangential multipliers are orders of magnitude smaller than the normal ones and would otherwise stall the iterations.
- **`print_convergence_criterion`** (local flag `PRINTING_OUTPUT`) selects between the plain `KRATOS_INFO` lines and the colored/bold output; when the `ProcessInfo` carries a `TABLE_UTILITY` (created by `AuxiliaryCreateConvergenceParameters` when `contact_settings.fancy_convergence_criterion` is `true`) each criterion adds its columns to the shared table in `Initialize` (guarded by `TABLE_IS_INITIALIZED`) and writes its values in `PostCriteria`.

### `BaseMortarConvergenceCriteria`

[Source](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_strategies/custom_convergencecriterias/base_mortar_criteria.h). Base of the five active-set criteria; its own `PreCriteria` and `PostCriteria` return `true`, they exist to keep the mortar quantities consistent with the current iterate. Defaults:

```json
{
    "name"                   : "base_mortar_criteria",
    "compute_dynamic_factor" : false,
    "gidio_debug"            : false,
    "pure_slip"              : false
}
```

The three booleans are stored in the local flags `COMPUTE_DYNAMIC_FACTOR`, `IO_DEBUG` and `PURE_SLIP`; the Python factory passes them positionally from `contact_settings.compute_dynamic_factor`, `contact_settings.gidio_debug` and the pure-slip detection (`"PureSlip" in mortar_type`, or `AuxiliaryPureSlipCheck` returning `true` when the sum of the `FRICTION_COEFFICIENT` of all properties is zero).

**`PreCriteria`** (lines 168-232), executed before the system is built in every iteration:

1. If `ProcessInfo[CONSIDER_NORMAL_VARIATION]` is not `NO_DERIVATIVES_COMPUTATION`, `ComputeNodesMeanNormalModelPartWithPairedNormal` recomputes the nodal normals of the `Contact` sub-model-part and refreshes the paired normal stored in every `PairedCondition`.
2. If the model part `Is(SLIP)` (frictional problem), the nodal tangents are updated by `MortarUtilities::ComputeNodesTangentModelPart`: from the tangential multiplier when the model part has `VECTOR_LAGRANGE_MULTIPLIER` and the problem is not pure slip, from `WEIGHTED_SLIP` otherwise. Tangents must be updated even when the normal is constant.
3. If `ProcessInfo[ADAPT_PENALTY]` is `true` or the problem is dynamic (`VELOCITY` is a nodal solution-step variable), the weighted gap is reset (`ResetWeightedGap`) and recomputed with `ContactUtilities::ComputeExplicitContributionConditions("ComputingContact")`.
4. In the dynamic case with `COMPUTE_DYNAMIC_FACTOR`, `ComputeDynamicFactorProcess` is executed on `Contact` (it computes the nodal `DYNAMIC_FACTOR` used to scale the gap).
5. With `ADAPT_PENALTY`, `AALMAdaptPenaltyValueProcess` updates the nodal `INITIAL_PENALTY` (adaptive augmented Lagrangian, thesis App. D.4.3.1).

**`PostCriteria`** (lines 235-290), executed after the update of the database:

1. Copies the current `WEIGHTED_GAP` of every `Contact` node into buffer position 1 (needed by the frictional formulation and by `MPCContactCriteria`).
2. Resets `WEIGHTED_GAP` (and `WEIGHTED_SLIP`) and recomputes them with `ComputeExplicitContributionConditions`, so that the active-set check that follows in the derived class is performed with the gap of the *updated* configuration.
3. If `IO_DEBUG` is set, a `GidIO` frame labelled with `NL_ITERATION_NUMBER` is written with the flags `INTERFACE`, `ACTIVE`, `SLAVE`, `ISOLATED` (and `SLIP`), the nodal `NORMAL`, `DYNAMIC_FACTOR`, `AUGMENTED_NORMAL_CONTACT_PRESSURE` (and `AUGMENTED_TANGENT_CONTACT_PRESSURE`), `DISPLACEMENT`, `VELOCITY` / `ACCELERATION` when present, the multipliers (`LAGRANGE_MULTIPLIER_CONTACT_PRESSURE` or `VECTOR_LAGRANGE_MULTIPLIER`) and `WEIGHTED_GAP` (and `WEIGHTED_SLIP`). This is the fastest way to debug an oscillating active set.

`Initialize` of the derived classes adds the columns `ACTIVE SET CONV` (and `SLIP/STICK CONV` for the frictional ones) to the shared table.

### The five active-set criteria

The derived criteria override `PostCriteria`, call the base version first and then call the corresponding function of the `ActiveSetUtilities` namespace ([`active_set_utilities.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_utilities/active_set_utilities.cpp)). Each function loops over the `SLAVE` nodes of the `Contact` sub-model-part, evaluates the nodal complementarity function, toggles the flags and returns the number of nodes that changed state; the criterion converges when that number is zero. The result is stored in `ProcessInfo[ACTIVE_SET_CONVERGED]` (and `SLIP_SET_CONVERGED`), and the `ACTIVE_SET_COMPUTED` flag prevents the ALM criteria from evaluating the set twice in the same iteration. The set is only re-evaluated when `rModelPart.Is(INTERACTION) || NL_ITERATION_NUMBER == 1`, which is the mechanism behind the two loop modes of the strategy. Defaults of all five (only `"name"` differs):

```json
{
    "name"                        : "alm_frictionless_mortar_criteria",
    "print_convergence_criterion" : false
}
```

With $$\varepsilon$$ the penalty (nodal `INITIAL_PENALTY` if present, otherwise `ProcessInfo[INITIAL_PENALTY]`), $$k$$ the scale factor `SCALE_FACTOR`, $$c_\tau$$ the `TANGENT_FACTOR`, $$\tilde{g}_n$$ the weighted gap `WEIGHTED_GAP`, $$\tilde{\mathbf{g}}_\tau$$ the weighted slip `WEIGHTED_SLIP` and $$\mu$$ the nodal `FRICTION_COEFFICIENT`, the checks are:

| Criterion | `ActiveSetUtilities` function | Augmented normal pressure $$\bar{\lambda}_n$$ (stored in `AUGMENTED_NORMAL_CONTACT_PRESSURE`) | Active if | Stick/slip |
|---|---|---|---|---|
| `ALMFrictionlessMortarConvergenceCriteria` | `ComputeALMFrictionlessActiveSet(rModelPart)` | $$k\lambda_n + \varepsilon\tilde{g}_n$$ with $$\lambda_n$$ = `LAGRANGE_MULTIPLIER_CONTACT_PRESSURE` | $$\bar{\lambda}_n \lt 0$$; on activation $$\lambda_n \leftarrow \bar{\lambda}_n/k$$ | – |
| `ALMFrictionlessComponentsMortarConvergenceCriteria` | `ComputeALMFrictionlessComponentsActiveSet(rModelPart)` | $$k\,\mathbf{n}\cdot\boldsymbol{\lambda} + \varepsilon\tilde{g}_n$$ with $$\boldsymbol{\lambda}$$ = `VECTOR_LAGRANGE_MULTIPLIER` | $$\bar{\lambda}_n \lt 0$$; on activation $$\boldsymbol{\lambda} \leftarrow \mathbf{n}\,\bar{\lambda}_n/k$$ | – |
| `ALMFrictionalMortarConvergenceCriteria` | `ComputeALMFrictionalActiveSet(rModelPart, PureSlip, EchoLevel)` | $$k\,\mathbf{n}\cdot\boldsymbol{\lambda} + \varepsilon\tilde{g}_n$$ | $$\bar{\lambda}_n \lt 0$$; on activation $$\boldsymbol{\lambda} \leftarrow \mathbf{n}\bar{\lambda}_n/k + \bar{\boldsymbol{\lambda}}_\tau/k$$ (tangent part omitted if $$\mu = 0$$) | see below |
| `PenaltyFrictionlessMortarConvergenceCriteria` | `ComputePenaltyFrictionlessActiveSet(rModelPart)` | $$\varepsilon\tilde{g}_n$$ | $$\bar{\lambda}_n \lt 0$$ | – |
| `PenaltyFrictionalMortarConvergenceCriteria` | `ComputePenaltyFrictionalActiveSet(rModelPart, PureSlip, EchoLevel)` | $$\varepsilon\tilde{g}_n$$ | $$\bar{\lambda}_n \lt 0$$ | see below |

Deactivated nodes reset `WEIGHTED_SLIP` and the `SLIP` flag. This is exactly the nodal NCP function of thesis eq. 4.44, $$\mathcal{C}_{\lambda_n} = k\lambda_n - \max(0, k\lambda_n + \varepsilon\tilde{g}_n)$$, evaluated with the sign convention "negative pressure means compression".

For the **frictional** checks the tangential multiplier of an active node is $$\boldsymbol{\lambda}_\tau = \boldsymbol{\lambda} - (\mathbf{n}\cdot\boldsymbol{\lambda})\mathbf{n}$$ and the augmented tangent pressure (stored in `AUGMENTED_TANGENT_CONTACT_PRESSURE`) is

<p align="center">$$\bar{\boldsymbol{\lambda}}_\tau = k\boldsymbol{\lambda}_\tau + c_\tau\,\varepsilon\,\tilde{\mathbf{g}}_\tau \quad\text{(stick node)},\qquad \bar{\boldsymbol{\lambda}}_\tau = k\boldsymbol{\lambda}_\tau + c_{slip}\,c_\tau\,\varepsilon\,\tilde{\mathbf{g}}_\tau \quad\text{(slip node)}$$</p>

where $$c_{slip}$$ is `ProcessInfo[SLIP_AUGMENTATION_COEFFICIENT]` (0 if absent). The node is set to **slip** when $$\Vert\bar{\boldsymbol{\lambda}}_\tau\Vert / (-\mu\bar{\lambda}_n) \gt \theta$$ with $$\theta = 1$$ for a stick node and $$\theta = 1 - $$ `SLIP_THRESHOLD` for a node that is already slipping (a small hysteresis that avoids chattering), and to **stick** otherwise; in the slip case `AUGMENTED_TANGENT_CONTACT_PRESSURE` is overwritten by the Coulomb limit $$-\mu\bar{\lambda}_n\,\boldsymbol{\lambda}_\tau/\Vert\boldsymbol{\lambda}_\tau\Vert$$ (thesis eqs. 4.75–4.77). With `PureSlip` every active node is forced to `SLIP` and a stick detection only produces a warning. The penalty frictional version uses $$\bar{\boldsymbol{\lambda}}_\tau = c_\tau\varepsilon\tilde{\mathbf{g}}_\tau$$ and the stick condition $$\Vert\bar{\boldsymbol{\lambda}}_\tau\Vert \le -\mu\bar{\lambda}_n$$. The frictional utilities return an `array_1d<std::size_t, 2>` with the number of active-set changes and the number of stick/slip changes; the ALM criterion reports both separately, the penalty one requires the sum to vanish.

### `MPCContactCriteria`

[Source](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_strategies/custom_convergencecriterias/mpc_contact_criteria.h). Defaults: `{"name" : "mpc_contact_criteria"}`. In the MPC route there are no multipliers, so the contact pressure has to be recovered from the **reactions**:

- `PreCriteria` (lines 146-175) resets the non-historical `CONTACT_FORCE` of the `Contact` nodes, saves `WEIGHTED_GAP` in buffer position 1, recomputes the weighted gap (`ComputeWeightedGap`) and zeroes `NODAL_AREA`.
- `PostCriteria` (lines 187-380), when `NL_ITERATION_NUMBER > 0`: recomputes the weighted gap and `NODAL_AREA`; for every pair `ContactSub<i>` it maps `REACTION` from `MasterSubModelPart<i>` to `SlaveSubModelPart<i>` with a `SimpleMortarMapperProcessWrapper` configured as `{"distance_threshold" : 1.0e24, "update_interface" : false, "origin_variable" : "REACTION", "mapping_coefficient" : -1.0e0}` (the threshold is replaced by `ProcessInfo[DISTANCE_THRESHOLD]` when available), so that the slave reaction accumulates the (sign-inverted) master reaction. Then, for every `SLAVE` node, the normal contact force is $$f_n = \mathbf{R}\cdot\mathbf{n}$$, the contact pressure is $$p_n = f_n/$$`NODAL_MAUX`, the nodal gap is $$g = \tilde{g}_n/$$`NODAL_AREA`, and the node is **active** if $$p_n \lt -(\texttt{REACTION\_CHECK\_STIFFNESS\_FACTOR}\cdot E)$$ or $$g \lt 0$$, where $$E$$ is the `YOUNG_MODULUS` of the first element properties (the threshold is 0 when no Young modulus is found; the factor defaults to $$10^{-12}$$ in the criterion and is set to `mpc_contact_settings.reaction_check_stiffness_factor`, default $$10^{-10}$$, by `MPCContactProcess`). Active nodes store `CONTACT_FORCE` $$= f_n\mathbf{n}/$$`NODAL_PAUX` and `NORMAL_CONTACT_STRESS` $$= p_n$$. In the frictional case (model part `Is(SLIP)`) the tangential pressure $$\Vert\mathbf{R} - f_n\mathbf{n}\Vert/$$`NODAL_MAUX` is compared with $$-\mu f_n$$ to toggle `SLIP` and to store `TANGENTIAL_CONTACT_STRESS`.
- Finally every condition of `ComputingContact` whose slave nodes are all inactive is set inactive together with its constraint (`CONSTRAINT_POINTER`, error if missing): this is the "tension check" that releases the constraint when the interface would work in traction. Convergence is reached when no node changed its active (or slip) state.

### `MortarAndConvergenceCriteria`

[Source](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_strategies/custom_convergencecriterias/mortar_and_criteria.h). A specialization of the core `And_Criteria` that combines the *mechanical* criterion (first argument) and the *mortar* active-set criterion (second argument). Defaults:

```json
{
    "name"                        : "mortar_and_criteria",
    "print_convergence_criterion" : false
}
```

The constructor takes an optional `ConditionNumberUtility::Pointer`. What it adds to the plain `And_Criteria`:

- `Initialize` adds the leading column `ITER` (and `COND.NUM.` when a condition-number utility is given) to the `TABLE_UTILITY` table, then lets the children add theirs; `InitializeSolutionStep` prints the table header.
- `PostCriteria` writes `NL_ITERATION_NUMBER` in the row, delegates to `And_Criteria::PostCriteria`, and, if a condition-number utility exists, copies the system matrix and prints `GetConditionNumber(A)` (computed with the power-iteration eigenvalue solvers `power_iteration_highest_eigenvalue_solver` / `power_iteration_eigenvalue_solver` built by the Python factory when `contact_settings.condn_convergence_criterion` is `true`; this is expensive and meant for diagnosis only). On convergence the table footer is printed.

This is the object the Python factory always returns for the `contact_*` criteria, which is why the console shows the characteristic table with the displacement, multiplier and active-set columns (thesis Figs. 4.18 and 4.23).

### Displacement, residual, mixed and frictional families

These nine criteria are the *mechanical* half of the `MortarAndConvergenceCriteria`. They differ in **what** is measured (increment `Dx`, residual `b`, or residual for displacements and increment for the multipliers), and in **how** the multiplier block is split. Each family has a version without multipliers (used with the penalty formulations, where the only unknowns are displacements), with one multiplier block (scalar ALM, components ALM, mesh tying) and with three multiplier blocks (frictional ALM: normal, stick-tangent, slip-tangent).

| Measured quantity | No multipliers (penalty) | One multiplier block | Normal + stick + slip blocks |
|---|---|---|---|
| Increment (`Dx`) | `DisplacementContactCriteria` | `DisplacementLagrangeMultiplierContactCriteria` | `DisplacementLagrangeMultiplierFrictionalContactCriteria` |
| Residual (`b`) | `DisplacementResidualContactCriteria` | `DisplacementLagrangeMultiplierResidualContactCriteria` | `DisplacementLagrangeMultiplierResidualFrictionalContactCriteria` |
| Mixed (residual for $$\mathbf{u}$$, increment for $$\boldsymbol{\lambda}$$) | – (the factory reuses `DisplacementResidualContactCriteria`) | `DisplacementLagrangeMultiplierMixedContactCriteria` | `DisplacementLagrangeMultiplierMixedFrictionalContactCriteria` |

The residual versions store the norm of the first iteration (`INITIAL_RESIDUAL_IS_SET`) and use it as reference for the ratio; they require the RHS to be up to date, which is why `SetActualizeRHSFlag(True)` is set on the composite. The mixed versions are usually the most robust choice for ALM: the multiplier residual is the weak gap, which is not a good convergence indicator while the active set is still changing.

`displacement_contact_criteria.h`:

```json
{
    "name"                            : "displacement_contact_criteria",
    "ensure_contact"                  : false,
    "print_convergence_criterion"     : false,
    "displacement_relative_tolerance" : 1.0e-4,
    "displacement_absolute_tolerance" : 1.0e-9,
    "rotation_relative_tolerance"     : 1.0e-4,
    "rotation_absolute_tolerance"     : 1.0e-9
}
```

`displacement_residual_contact_criteria.h`:

```json
{
    "name"                                 : "displacement_residual_contact_criteria",
    "ensure_contact"                       : false,
    "print_convergence_criterion"          : false,
    "residual_relative_tolerance"          : 1.0e-4,
    "residual_absolute_tolerance"          : 1.0e-9,
    "rotation_residual_relative_tolerance" : 1.0e-4,
    "rotation_residual_absolute_tolerance" : 1.0e-9
}
```

`displacement_lagrangemultiplier_contact_criteria.h`:

```json
{
    "name"                                    : "displacement_lagrangemultiplier_contact_criteria",
    "ensure_contact"                          : false,
    "print_convergence_criterion"             : false,
    "displacement_relative_tolerance"         : 1.0e-4,
    "displacement_absolute_tolerance"         : 1.0e-9,
    "rotation_relative_tolerance"             : 1.0e-4,
    "rotation_absolute_tolerance"             : 1.0e-9,
    "contact_displacement_relative_tolerance" : 1.0e-4,
    "contact_displacement_absolute_tolerance" : 1.0e-9
}
```

`displacement_lagrangemultiplier_residual_contact_criteria.h`:

```json
{
    "name"                                 : "displacement_lagrangemultiplier_residual_contact_criteria",
    "ensure_contact"                       : false,
    "print_convergence_criterion"          : false,
    "residual_relative_tolerance"          : 1.0e-4,
    "residual_absolute_tolerance"          : 1.0e-9,
    "rotation_residual_relative_tolerance" : 1.0e-4,
    "rotation_residual_absolute_tolerance" : 1.0e-9,
    "contact_residual_relative_tolerance"  : 1.0e-4,
    "contact_residual_absolute_tolerance"  : 1.0e-9
}
```

`displacement_lagrangemultiplier_mixed_contact_criteria.h`:

```json
{
    "name"                                    : "displacement_lagrange_multiplier_mixed_contact_criteria",
    "ensure_contact"                          : false,
    "print_convergence_criterion"             : false,
    "residual_relative_tolerance"             : 1.0e-4,
    "residual_absolute_tolerance"             : 1.0e-9,
    "rotation_residual_relative_tolerance"    : 1.0e-4,
    "rotation_residual_absolute_tolerance"    : 1.0e-9,
    "contact_displacement_relative_tolerance" : 1.0e-4,
    "contact_displacement_absolute_tolerance" : 1.0e-9
}
```

`displacement_lagrangemultiplier_frictional_contact_criteria.h`:

```json
{
    "name"                                                     : "displacement_lagrangemultiplier_frictional_contact_criteria",
    "ensure_contact"                                           : false,
    "pure_slip"                                                : false,
    "print_convergence_criterion"                              : false,
    "displacement_relative_tolerance"                          : 1.0e-4,
    "displacement_absolute_tolerance"                          : 1.0e-9,
    "rotation_relative_tolerance"                              : 1.0e-4,
    "rotation_absolute_tolerance"                              : 1.0e-9,
    "contact_displacement_relative_tolerance"                  : 1.0e-4,
    "contact_displacement_absolute_tolerance"                  : 1.0e-9,
    "frictional_stick_contact_displacement_relative_tolerance" : 1.0e-4,
    "frictional_stick_contact_displacement_absolute_tolerance" : 1.0e-9,
    "frictional_slip_contact_displacement_relative_tolerance"  : 1.0e-4,
    "frictional_slip_contact_displacement_absolute_tolerance"  : 1.0e-9,
    "ratio_normal_tangent_threshold"                           : 1.0e-4
}
```

`displacement_lagrangemultiplier_residual_frictional_contact_criteria.h`:

```json
{
    "name"                                                 : "displacement_lagrangemultiplier_ressidual_frictional_contact_criteria",
    "ensure_contact"                                       : false,
    "pure_slip"                                            : false,
    "print_convergence_criterion"                          : false,
    "residual_relative_tolerance"                          : 1.0e-4,
    "residual_absolute_tolerance"                          : 1.0e-9,
    "rotation_residual_relative_tolerance"                 : 1.0e-4,
    "rotation_residual_absolute_tolerance"                 : 1.0e-9,
    "contact_residual_relative_tolerance"                  : 1.0e-4,
    "contact_residual_absolute_tolerance"                  : 1.0e-9,
    "frictional_stick_contact_residual_relative_tolerance" : 1.0e-4,
    "frictional_stick_contact_residual_absolute_tolerance" : 1.0e-9,
    "frictional_slip_contact_residual_relative_tolerance"  : 1.0e-4,
    "frictional_slip_contact_residual_absolute_tolerance"  : 1.0e-9
}
```

`displacement_lagrangemultiplier_mixed_frictional_contact_criteria.h`:

```json
{
    "name"                                                     : "displacement_lagrangemultiplier_mixed_frictional_contact_criteria",
    "ensure_contact"                                           : false,
    "pure_slip"                                                : false,
    "print_convergence_criterion"                              : false,
    "residual_relative_tolerance"                              : 1.0e-4,
    "residual_absolute_tolerance"                              : 1.0e-9,
    "rotation_residual_relative_tolerance"                     : 1.0e-4,
    "rotation_residual_absolute_tolerance"                     : 1.0e-9,
    "contact_displacement_relative_tolerance"                  : 1.0e-4,
    "contact_displacement_absolute_tolerance"                  : 1.0e-9,
    "frictional_stick_contact_displacement_relative_tolerance" : 1.0e-4,
    "frictional_stick_contact_residual_relative_tolerance"     : 1.0e-9,
    "frictional_slip_contact_displacement_relative_tolerance"  : 1.0e-4,
    "frictional_slip_contact_residual_relative_tolerance"      : 1.0e-9,
    "ratio_normal_tangent_threshold"                           : 1.0e-4
}
```

Note: in the mixed frictional block the keys `frictional_stick_contact_residual_relative_tolerance` and `frictional_slip_contact_residual_relative_tolerance` play the role of the *absolute* tolerances of the tangential blocks (their default is $$10^{-9}$$); the Python factory maps them accordingly (`FSTCR_AT`, `FSLCR_AT`). In the frictional criteria the convergence of the multiplier block is

<p align="center">$$\text{conv}_\lambda = \text{conv}(\lambda_n) \;\wedge\; \big[\text{conv}(\lambda_{\tau,stick}) \vee r_{stick}/r_n \le \theta_{nt}\big] \;\wedge\; \big[\text{conv}(\lambda_{\tau,slip}) \vee r_{slip}/r_n \le \theta_{nt}\big]$$</p>

with $$\theta_{nt}$$ the `ratio_normal_tangent_threshold` and $$r$$ the absolute norms (thesis eq. 4.78 checks the four residuals $$\mathbf{r}_u$$, $$\mathbf{r}_{\lambda_n}$$, $$\mathbf{r}_{\lambda_\tau^{sl}}$$, $$\mathbf{r}_{\lambda_\tau^{st}}$$ separately for the same reason).

### `ContactErrorMeshCriteria`

[Source](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_strategies/custom_convergencecriterias/contact_error_mesh_criteria.h). Used only by the adaptive-remeshing analysis stage (`adaptive_remeshing/adaptative_remeshing_contact_structural_mechanics_utilities.py`, when the error criterion is `adaptative_remesh_criteria`), never by the plain solvers; it needs the `MeshingApplication`. Defaults:

```json
{
    "name"                 : "contact_error_mesh_criteria",
    "error_mesh_tolerance" : 5.0e-3,
    "error_mesh_constant"  : 5.0e-3,
    "compute_error_extra_parameters":
    {
        "penalty_normal"       : 1.0e4,
        "penalty_tangential"   : 1.0e4,
        "echo_level"           : 0
    }
}
```

`PostCriteria` flags the nodes and conditions of `Contact` with `CONTACT`, runs `ContactSPRErrorProcess<2>` or `<3>` (a superconvergent-patch-recovery error estimator that accounts for the contact tractions through the two penalties) and compares `ProcessInfo[ERROR_RATIO]` with `error_mesh_tolerance`. Returning `true` means that the discretization error is below the tolerance and **no remeshing is required**. See [Adaptive remeshing](../Examples/Adaptive_Remeshing.html).

### `MeshTyingMortarConvergenceCriteria`

[Source](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_strategies/custom_convergencecriterias/mesh_tying_mortar_criteria.h). Defaults: `{"name" : "mesh_tying_mortar_criteria"}`. Mesh tying has no inequality, hence no active set; the class only participates in the table printing and always returns `true`. It is what `GetMortarCriteria` returns when `"MeshTying" in mortar_type`, so that the same `MortarAndConvergenceCriteria` composition works for tying and for contact (see [Mesh tying](../Theory/Mesh_Tying.html)).

## How the Python factory composes the criteria

`ContactConvergenceCriteriaFactory` ([`contact_convergence_criteria_factory.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/python_scripts/contact_convergence_criteria_factory.py)) receives the parameters assembled by `auxiliary_methods_solvers.AuxiliaryCreateConvergenceParameters`: the standard keys of `solver_settings` (`convergence_criterion`, `rotation_dofs`, `echo_level`, `displacement_*` and `residual_*` tolerances) plus the contact-specific keys copied from `contact_settings` (all the `rotation_*`, `contact_*`, `frictional_*` tolerances, `ratio_normal_tangent_threshold`, `mortar_type`, `condn_convergence_criterion`, `print_convergence_criterion`, `ensure_contact`, `frictional_decomposed`, `compute_dynamic_factor`, `gidio_debug`). The composition depends on the `convergence_criterion` string:

| `convergence_criterion` | `mortar_type` contains `Penalty` | `ALMContactFrictional*` and `frictional_decomposed` | Otherwise (scalar / components ALM, mesh tying) |
|---|---|---|---|
| `contact_displacement_criterion` | `DisplacementContactCriteria` | `DisplacementLagrangeMultiplierFrictionalContactCriteria` | `DisplacementLagrangeMultiplierContactCriteria` |
| `contact_residual_criterion` | `DisplacementResidualContactCriteria` | `DisplacementLagrangeMultiplierResidualFrictionalContactCriteria` | `DisplacementLagrangeMultiplierResidualContactCriteria` |
| `contact_mixed_criterion` | `DisplacementResidualContactCriteria` (no rotation tolerances) | `DisplacementLagrangeMultiplierMixedFrictionalContactCriteria` | `DisplacementLagrangeMultiplierMixedContactCriteria` |
| `contact_and_criterion` | `KM.AndCriteria(DisplacementResidualContactCriteria, DisplacementContactCriteria)` | (same as "otherwise") | `KM.AndCriteria(DisplacementLagrangeMultiplierResidualContactCriteria, DisplacementLagrangeMultiplierContactCriteria)` |
| `contact_or_criterion` | `KM.OrCriteria(...)` of the same pair | (same as "otherwise") | `KM.OrCriteria(...)` of the same pair |
| `adaptative_remesh_criteria` | `None` (the adaptive-remeshing analysis builds `ContactErrorMeshCriteria` itself) | | |
| any other value (`displacement_criterion`, `residual_criterion`, `and_criterion`, `or_criterion`, …) | The StructuralMechanicsApplication `convergence_criteria_factory` builds the criterion; if `mortar_type` contains `ALMContact` or `MeshTying` it is combined with the mortar criterion in a plain `KM.AndCriteria` (no table) with `SetActualizeRHSFlag(True)` | | |

For every `contact_*` value the mechanical criterion is then wrapped as `CSMA.MortarAndConvergenceCriteria(mechanical, Mortar, print_convergence_criterion, condition_number_utility)`, where `Mortar = GetMortarCriteria()` is chosen from `mortar_type`:

| `mortar_type` | Mortar criterion |
|---|---|
| `ALMContactFrictionless` | `ALMFrictionlessMortarConvergenceCriteria(print, compute_dynamic_factor, gidio_debug)` |
| `ALMContactFrictionlessComponents` | `ALMFrictionlessComponentsMortarConvergenceCriteria(...)` |
| contains `ALMContactFrictional` | `ALMFrictionalMortarConvergenceCriteria(pure_slip, ...)` |
| `PenaltyContactFrictionless` | `PenaltyFrictionlessMortarConvergenceCriteria(...)` |
| contains `PenaltyContactFrictional` | `PenaltyFrictionalMortarConvergenceCriteria(pure_slip, ...)` |
| contains `MeshTying` | `MeshTyingMortarConvergenceCriteria()` |

`pure_slip` is `True` when `mortar_type` contains `PureSlip`, otherwise it is detected by `AuxiliaryPureSlipCheck` (all `FRICTION_COEFFICIENT` values of the properties are zero). The MPC solvers do not use this factory: `MPCContactStaticSolver._CreateConvergenceCriterion` builds the structural criterion and combines it with `MPCContactCriteria()` in a `KM.AndCriteria`.

## Relation to the thesis algorithms

The code above is the implementation of thesis Algorithm 2 (frictionless, §4.3.3.5) and Algorithm 3 (frictional, §4.3.4.4), whose *active-set strategy* is the primal-dual active-set / semi-smooth Newton method of §4.3.3.6 and §4.3.4.5. In pseudo-code, annotated with the classes that perform each line:

```
Algorithm 2 (thesis) — frictionless contact               Implementation
1  t = 0, i = 0; u^0 = 0, lambda^0 = 0                     StructuralMechanicsAnalysis, AuxiliaryAddDofs
2  Initialize the active set A_1^0, I_1^0                  search (ACTIVE_CHECK_FACTOR) in SearchBaseProcess
3  while t < t_end:
4     t = t + dt, i = i + 1; Delta u_1^i = 0                strategy InitializeSolutionStep / Predict
5     search contact pairs, update pairs and active set    SearchBaseProcess.ExecuteInitializeSolutionStep
6     conv = false
7     while not conv:                                      ResidualBasedNewtonRaphsonContactStrategy::BaseSolveSolutionStep
8        PreCriteria: normals, tangents, gap, penalty      BaseMortarConvergenceCriteria::PreCriteria
9        solve the linearized system (thesis 4.3.3.4.3)    BuilderAndSolver::BuildAndSolve (+ MixedULMLinearSolver)
10       u^i_{n+1} = u^i_n + Delta u, lambda likewise      UpdateDatabase
11       update active set with threshold (eq. 4.41-4.42)  ActiveSetUtilities::Compute*ActiveSet via PostCriteria
              A := {j : k*lambda_n + eps*g_n < 0}
              I := {j : k*lambda_n + eps*g_n >= 0}
12       check ||r_u|| < tol_u, ||r_lambda|| < tol_lambda   Displacement*ContactCriteria::PostCriteria
13       conv = (A, I unchanged) and residuals converged   MortarAndConvergenceCriteria (And)
14    FinalizeSolutionStep; output                          strategy, SearchBaseProcess.ExecuteFinalizeSolutionStep
```

```
Algorithm 3 (thesis) — frictional contact                  Differences with respect to Algorithm 2
2  also initialize the slip/stick sets A_sl, A_st          SLIP flag (search pre-activation, slip_step_reset_frequency)
8  tangents are recomputed every iteration                 MortarUtilities::ComputeNodesTangentModelPart in PreCriteria
11 update A, I (eq. 4.73-4.74) and then A_sl, A_st         ActiveSetUtilities::ComputeALMFrictionalActiveSet
      with the frictional threshold F = mu*(k n.lambda + eps*g_n)   (eq. 4.76) and the
      tangent stress t = ||k tau.lambda + eps_tau u_tau||           (eq. 4.77)
12 four residuals: r_u, r_lambda_n, r_lambda_tau^sl, r_lambda_tau^st  (eq. 4.78)
                                                            DisplacementLagrangeMultiplier*FrictionalContactCriteria
13 conv = (A, I, A_sl, A_st unchanged) and residuals        ALMFrictionalMortarConvergenceCriteria reports
                                                            ACTIVE SET CONV and SLIP/STICK CONV separately
```

The only structural deviation from the thesis algorithms is the `simplified_semi_smooth_newton` mode, in which line 11 is executed after a complete Newton loop instead of inside every iteration; it trades robustness (the active set cannot oscillate within the Newton iterations) for cost (each active-set update costs a full Newton solve) and is the mode of choice when the full semi-smooth Newton chatters, see [Tips, troubleshooting and limitations](../Usage/Tips_Troubleshooting_And_Limitations.html).

<p align="center"><img src="images/csma_active_set_flowchart.svg" alt="Flowchart of the semi-smooth Newton active-set loop: PreCriteria, build and solve, update, PostCriteria with ActiveSetUtilities, INTERACTION branch and inner loop" width="1000"/></p>
<p align="center"><em>Figure: the active-set flowchart, mapping thesis Algorithms 2 and 3 to the strategy, the criteria and <code>ActiveSetUtilities</code>.</em></p>

The theory behind the thresholds (augmented pressure, NCP functions, Figs. 4.19 and 4.24 of the thesis) is developed in [Frictionless contact](../Theory/Frictionless_Contact.html) and [Frictional contact](../Theory/Frictional_Contact.html); the linear systems solved in line 9 and the way the multipliers are condensed are in [Builder and solvers and linear solvers](Builder_And_Solvers_And_Linear_Solvers.html).

## Quick reference: which settings reach which class

| `solver_settings` key | Consumed by |
|---|---|
| `solving_strategy_settings.type` (`newton_raphson`, `line_search`, `arc_length`) | `ContactStaticMechanicalSolver._CreateSolutionStrategy` → `AuxiliaryNewton` / `AuxiliaryLineSearch` |
| `max_iteration`, `compute_reactions`, `reform_dofs_at_each_step`, `move_mesh_flag` | Strategy constructors (positional arguments) |
| `convergence_criterion`, `displacement_*`, `residual_*`, `rotation_dofs`, `echo_level` | `ContactConvergenceCriteriaFactory` |
| `contact_settings.simplified_semi_smooth_newton` | `INTERACTION` flag of the computing model part (contact solvers) or of its `ProcessInfo` (MPC solvers) |
| `contact_settings.adaptative_strategy`, `split_factor`, `max_number_splits`, `inner_loop_iterations` | `ResidualBasedNewtonRaphsonContactStrategy` |
| `contact_settings.fancy_convergence_criterion` | Creates the `TABLE_UTILITY` (`TableStreamUtility`) read by every criterion |
| `contact_settings.condn_convergence_criterion` | `ConditionNumberUtility` passed to `MortarAndConvergenceCriteria` |
| `contact_settings.print_convergence_criterion`, `ensure_contact`, `frictional_decomposed`, `compute_dynamic_factor`, `gidio_debug`, `contact_*` / `frictional_*` tolerances, `ratio_normal_tangent_threshold` | The mechanical and mortar criteria through the factory |
| `contact_settings.silent_strategy` | `ContactStaticMechanicalSolver.Initialize` sets the strategy echo level to 0 right after the base `Initialize` has created it |
| `mpc_contact_settings.inner_loop_iterations`, `update_each_nl_iteration`, `enforce_ntn`, `simplified_semi_smooth_newton` | `ResidualBasedNewtonRaphsonMPCContactStrategy` |

The complete list of keys with their defaults is in the [Solver settings reference](../Usage/Solver_Settings_Reference.html).
