---
title: Tips Troubleshooting And Limitations
keywords: troubleshooting, convergence, penalty, scale factor, active set chattering, normals, master slave, limitations, known issues
tags: [usage, troubleshooting, tips, limitations, known issues, FAQ]
sidebar: contact_structural_mechanics_application
summary: Practical guidance for setting up and debugging contact simulations — parameter calibration, convergence problems, search problems, model set-up pitfalls, performance — followed by the known limitations and defects of the current implementation.
---

> **Sources.** Thesis §4.3.3.3 (ALM parameter calibration, Tables 4.1–4.2, Figs. 4.13–4.14), §4.3.3.6 and §4.3.4.5 (active-set strategy), §4.4 (search), §8.3 and §8.7 (conclusions and remarks); the code of the application as documented in the [Implementation](../Implementation/Architecture.html) pages. The defects listed at the end were found while writing this documentation and are reported here without modifying the code.

## Choosing the formulation

| Need | Use | Why |
|---|---|---|
| Robust default for frictionless problems | `alm_contact_process`, `contact_type: Frictionless` (scalar multiplier) | Exact constraint satisfaction, smallest extra system, no parameter tuning required in most cases. |
| Iterative linear solvers / many contact DoFs | `FrictionlessComponents` with the default `use_mixed_ulm_solver` | The dual multipliers are condensed and the displacement block can be solved with AMGCL (block size = dimension). |
| Friction | `Frictional` (or `FrictionalPureSlip` when stick states do not matter) | Coulomb law with consistent slip linearisation; `PureSlip` avoids stick/slip chattering. |
| Impact, explicit dynamics | `explicit_penalty_contact_process` | No Newton loop; the penalty limits the penetration, `delta_time_factor_for_contact` keeps the step stable. |
| Displacement-only system, no multipliers, legacy MPC solvers | `mpc_contact_process` | Constraints instead of multipliers; less accurate on non-matching meshes (thesis Fig. D.5). |
| Gluing non-matching meshes | `mesh_tying_process` | Equality constraint, no active set. |

Details in [Contact process settings](Contact_Process_Settings_Reference.html).

## Convergence troubleshooting

**The active set never settles (chattering).** Rows alternate `Achieved` / `Not achieved` in the `ACTIVE SET CONV` column.
- Reduce the load increment (`time_step`) or enable `adaptative_strategy` so that the step is split automatically (`split_factor`, `max_number_splits`).
- Check the penalty/scale factor: too large an $$\varepsilon$$ makes the augmented pressure $$\bar\lambda_n = k\lambda_n + \varepsilon\tilde g_n$$ over-sensitive to tiny gaps. The thesis experiment (Tables 4.1–4.2, Fig. 4.14) shows that $$\varepsilon$$ always worsens the condition number while $$k$$ has an optimum near $$\varepsilon \approx k \approx 10\,E/h$$; start from the automatic values (`stiffness_factor` 1–10) and change `penalty_scale_factor` only if needed.
- Try `simplified_semi_smooth_newton: true` (the sets are frozen during an inner Newton loop) with `inner_loop_iterations` 5–10, or the opposite if you were using it.
- For frictional problems increase `slip_threshold` (hysteresis of the stick/slip decision) or use `FrictionalPureSlip`.
- `predict_correct_lagrange_multiplier` (advanced search) gives a better initial multiplier after the search.

**Residual ratio stalls or diverges.**
- Large penetrations at the first iteration: reduce `active_check_factor`/`search_factor` only if wrong pairs are created, otherwise reduce the step.
- Rotating or strongly deforming interfaces: use `normal_variation: "nodal_elemental_derivatives"` (the `NV` conditions include the normal derivatives in the tangent, thesis §4.6) or at least `no_derivatives_computation_with_normal_update`.
- Very stiff bodies against soft ones: the automatic $$\varepsilon = E_{mean}/h_{mean}$$ uses the mean Young modulus of the interface; with a rigid body set `manual_ALM` and choose $$\varepsilon$$ from the *soft* body.
- Use `line_search` (`solving_strategy_settings.type`) — the contact line search fits separate parabolas for the displacement and multiplier increments.

**The multiplier residual is large while the displacement converges.** Loosen `contact_residual_relative_tolerance` slightly or switch to `contact_mixed_criterion` (displacement residual + multiplier increment), which is less sensitive to the scaling of the multipliers.

**Singular or ill-conditioned matrix.**
- Slave nodes whose pairs were all deactivated carry multipliers without equations: the block builder (`builder_and_solver_settings.type: "block"`) fixes them (`ISOLATED` flag); the elimination builder does not — use `block`.
- Multipliers on Dirichlet nodes: the elimination builder fixes the multiplier of a slave node whose displacement is fixed; avoid slave nodes with fully prescribed displacement when using the block builder.
- Check the condition number with `condn_convergence_criterion` (thesis §4.3.3.3) and tune $$k$$.

**Nothing comes into contact / everything is active from the start.** Inspect the flags `ACTIVE`, `SLAVE`, `MASTER` and `WEIGHTED_GAP` in the output of the first step; `ensure_contact: true` turns "no active node" into an error; `debug_mode` prints the total contact force against the applied load.

## Search troubleshooting

- **Inverted normals.** `NormalCheckProcess` detects and corrects conditions whose normal points inwards (`normal_check_proportion`); if a body still does not see the other one, check the orientation of the skin conditions in the mesh (the tests `test_check_normals_process.py` show the expected behaviour).
- **Master or slave?** The slave side is integrated (mortar operators, multipliers). Prefer as slave the surface that is finer, more curved or expected to receive the higher pressure; a coarse rigid tool is a natural master. Only nodes can carry multipliers, so a slave surface with very few nodes gives a poor pressure resolution.
- **Pairs not found.** Increase `search_factor` (search radius as a multiple of the condition size) or switch to `in_box_with_obb` / `octree_with_obb`; `adapt_search` rescales the factors with the relative mesh size of both sides; for fast-moving bodies use `dynamic_search`.
- **Too many pairs / far pairs activated.** Enable `consider_gap_threshold`, reduce `bounding_box_factor`, or lower `active_check_factor`; `normal_orientation_threshold` removes pairs whose normals are almost parallel (same-facing facets).
- **Wrong gap on curved surfaces.** Keep `check_gap: check_mapping` (the mortar-mapped consistent gap, thesis Algorithm 8) rather than `direct_check`.
- **Pairs frozen.** `database_step_update` controls how often the search runs; for large sliding keep it at 1.
- **Self-contact.** Leave `assume_master_slave` empty; conditions sharing master and slave nodes are deactivated by design ([Self contact](../Contact_Search/Self_Contact.html)).

## Model set-up pitfalls

- The contact sub-model-parts must contain **conditions** (the skin); `InterfacePreprocessCondition` can create them from elements when `contact_model_part` points to a body, but explicit skin conditions are safer.
- Sub-model-part names in the JSON must match the mesh (`Structure.<name>`); the process creates `Contact`, `ContactSub<key>`, `MasterSubModelPart<key>`, `SlaveSubModelPart<key>` and `ComputingContact` — do not use these names yourself.
- Mixed 3D interfaces (triangles against quadrilaterals) are supported (`3D3N4N`, `3D4N3N`); higher-order geometries are not.
- Axisymmetric conditions exist only in 2D and need `THICKNESS` in the properties; they do not exist for the components formulation.
- Frictional problems need `buffer_size` ≥ 3 (set automatically) and `FRICTION_COEFFICIENT` per interface (`friction_coefficients`).
- `clear_storage` and `reform_dofs_at_each_step` are forced to `true`; do not rely on keeping the system matrix between steps.
- Restarts: the normal check is skipped when `IS_RESTARTED` is set; the pairs are rebuilt at the first step.
- Adaptive remeshing (`adaptative_remesh_criteria`, `contact_remesh_mmg_process.py`) requires the `MeshingApplication` compiled with MMG; without it the corresponding solver raises an error.

## Performance

- Assembly: the generated frictional conditions are large; the frictional ALM `.cpp` alone has ~170 000 lines and dominates the compile time (excluded from unity builds). At run time the cost per pair is the exact segmentation plus the derivative computation; `integration_order: 2` is enough for linear geometries.
- Linear solver: for vector multipliers keep `use_mixed_ulm_solver` and an AMGCL inner solver; for the scalar formulation a direct solver on the mixed (saddle point) system is the safest choice; AMGCL on the mixed system needs the block builder and may struggle with the zero diagonal block.
- Search: `in_radius_with_obb` is the default trade-off; `octree_with_obb` scales better for very large surfaces; increase `bucket_size` for large KD-trees.
- Output: `clear_inactive_for_post` removes inactive pairs before writing; `gidio_debug` and `debug_mode` write a file per iteration — disable them in production.
- Parallelism: OpenMP only (the application has no MPI implementation).

## Known limitations

- No MPI / distributed memory support; no higher-order (quadratic) contact geometries; 2D axisymmetry only for the scalar and frictional ALM and the penalty families.
- Frictional laws: only Coulomb is wired into the conditions; the `frictional_law` key and the `TrescaFrictionalLaw` class exist but are not used by the conditions (`FRICTIONAL_LAW` variable registered, never read).
- Contact with beams/shells relies on their skin conditions and displacement DoFs only; rotations are not coupled.
- The mesh-tying `3D3N4N` / `3D4N3N` prototypes are registered with matching geometries (the created pairs are correct).
- `enforce_ntn` of the MPC strategy is accepted but its implementation is commented out.
- The contact remeshing examples formerly published in the Examples repository (`mmg_remeshing_examples/use_cases/contact_*`) are no longer available; the remeshing pages of this documentation rely on the thesis and the code.

## Known defects (documentation notes, code unchanged)

| Where | Issue | Effect |
|---|---|---|
| `custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_residual_frictional_contact_criteria.h` | Registered name `displacement_lagrangemultiplier_ressidual_frictional_contact_criteria` (double *s*). | Only relevant when creating the criterion by name through the registry. |
| `python_scripts/replace_properties_process.py` | The default-settings JSON string contains `"reinitialize_entities" : false.` (period instead of comma). | The process fails to parse its defaults; it is not used by the contact processes. |
| `custom_conditions/ALM_frictional_mortar_contact_axisym_condition.h`, `penalty_frictionless_mortar_contact_axisym_condition.h` | `MatrixSize` is re-declared with a value different from the base family. | Harmless: the resize helpers use the base constant. |
| `custom_python/contact_structural_mechanics_python_application.cpp` | `ACTIVE_CHECK_FACTOR` registered twice; `TANGENT_FACTOR` (a core variable) re-registered; `CONSTRAINT_POINTER`, `PARENT_ELEMENT`, `FRICTIONAL_LAW` not exposed. | Cosmetic; the pointer variables are not usable from Python. |
| `automatic_differentiation/penalty_frictionless_mortar_condition/README.md` | Names the frictionless-ALM script and output instead of the penalty ones. | Documentation only; the correct files are listed in [Automatic differentiation](../Theory/Automatic_Differentiation.html). |
| `automatic_differentiation/ALM_frictional_mortar_condition/alm_frictional_mortar_contact_condition.tex` | Unfilled template stub; the components `.tex` is byte-identical to the frictionless one. | Documentation only. |
| `tests/test_ContactStructuralMechanicsApplication.py` | `TALMHertzSphereTestContact` / `TComponentsALMHertzSphereTestContact` (axisymmetric, memory), `TALMIroningTestContact`, `TALMIroningDieTestContact`, `TMultiLayerContactTest` are commented out; one Components elimination test is skipped on Windows. | Reduced coverage of the validation suite. |
| `automatic_differentiation/*` | The generators need sympy 1.2 (symbol-to-function conversion removed in 1.3). | Regeneration requires a dedicated environment. |
| `README.md` of the application | Advertised test counts are outdated (the runner registers 122 active tests, 91 gtests). | Documentation only. |

Please report new issues on the [Kratos issue tracker](https://github.com/KratosMultiphysics/Kratos/issues) with the label of the application.
