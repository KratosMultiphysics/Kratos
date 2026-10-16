---
title: Frictional Laws and MPC Constraint
keywords: contact, friction, Coulomb, Tresca, frictional law, FRICTION_COEFFICIENT, multipoint constraint, MPC, ContactMasterSlaveConstraint, CONSTRAINT_POINTER, master-slave elimination, Kratos
tags: [contact, implementation, friction, mpc, constraints, ContactStructuralMechanicsApplication]
sidebar: contact_structural_mechanics_application
summary: The FrictionalLaw hierarchy (Coulomb and Tresca thresholds and their linearization), the current status of its integration in the conditions, and the constraint-based (MPC) contact route built on ContactMasterSlaveConstraint, MPCMortarContactCondition, the MPC strategy and the MPC criterion.
---

> **Sources.** Thesis §4.2.5 (frictional models, p. 96, Fig. 4.10), §4.3.4.1.1 (Coulomb's law, p. 113, eqs. 4.45), App. D.5 (multipoint constraints, master–slave elimination, pp. 320–321); code: `custom_frictional_laws/frictional_law.{h,cpp}`, `frictional_law_with_derivative.{h,cpp}`, `coulomb_frictional_law.{h,cpp}`, `tresca_frictional_law.{h,cpp}`, `custom_python/add_custom_frictional_laws_to_python.cpp`, `custom_master_slave_constraints/contact_master_slave_constraint.{h,cpp}`, `custom_conditions/mpc_mortar_contact_condition.{h,cpp}`, `custom_processes/mpc_contact_search_process.cpp`, `custom_strategies/custom_strategies/residualbased_newton_raphson_mpc_contact_strategy.h`, `custom_strategies/custom_convergencecriterias/mpc_contact_criteria.h`, `python_scripts/mpc_contact_process.py`, `python_scripts/mpc_contact_structural_mechanics_static_solver.py`, `custom_processes/alm_fast_init_process.cpp`.

This page groups two small but self-contained parts of the application. The first is the **frictional law** hierarchy of `custom_frictional_laws/`, a set of classes that encapsulate the tangential threshold $$\mathscr{F}$$ of a friction model (Coulomb, Tresca) and its linearization, intended as the plug-in point for other friction models in the mortar conditions. The second is the **multipoint-constraint (MPC) contact route**: an alternative to Lagrange multipliers and penalties in which the mortar coupling between slave and master displacements is imposed as a Kratos `MasterSlaveConstraint` (`ContactMasterSlaveConstraint`, `custom_master_slave_constraints/`) that the contact search creates, the `MPCMortarContactCondition` updates and the builder-and-solver eliminates. Both are documented here because both are *interfaces to the core machinery* rather than formulations of their own: the mathematics of the frictional contact conditions is in [Frictional contact](../Theory/Frictional_Contact.html), the master–slave elimination in [Constrained optimisation methods](../Theory/Constrained_Optimisation_Methods.html).

## Frictional laws

### Hierarchy and interface

```
FrictionalLaw                                                     frictional_law.{h,cpp}              (not templated)
└── FrictionalLawWithDerivative<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>
                                                                  frictional_law_with_derivative.{h,cpp}
    ├── CoulombFrictionalLaw<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>
    │                                                             coulomb_frictional_law.{h,cpp}
    └── TrescaFrictionalLaw<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>
                                                                  tresca_frictional_law.{h,cpp}
```

All classes are exported with `KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION)` and compiled into the core library (`custom_frictional_laws/*.cpp` is one of the globs of `CMakeLists.txt`). The base class ([`frictional_law.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_frictional_laws/frictional_law.h)) is documented as "the base class for frictional laws … this class does nothing, define derived frictional laws in order to make use of it". Its interface:

| Method | Class | Behavior |
|---|---|---|
| `double GetFrictionCoefficient(const Node& rNode, const PairedCondition& rCondition, const ProcessInfo& rCurrentProcessInfo)` | `FrictionalLaw` (virtual) | Looks up `FRICTION_COEFFICIENT` in this order: `rCondition.GetProperties()`, the nodal non-historical data of `rNode`, the `ProcessInfo`; returns 0.0 if none has it |
| `double GetThresholdValue(const Node&, const PairedCondition&, const ProcessInfo&)` | `FrictionalLaw` (virtual) | The tangential threshold $$\mathscr{F}$$ of the node. The base implementation raises `KRATOS_ERROR` ("You are calling to the base class method GetThresholdValue, check your frictional law declaration") |
| `double GetDerivativeThresholdValue(const Node&, const PairedCondition&, const ProcessInfo&, const DerivativeDataType& rDerivativeData, const MortarConditionMatrices& rMortarConditionMatrices, IndexType IndexDerivative, IndexType IndexNode)` | `FrictionalLawWithDerivative` (virtual) | Directional derivative of the threshold of slave node `IndexNode` with respect to the DoF `IndexDerivative` of the pair condition (displacements of the slave and master nodes first, then the multipliers). The base implementation raises `KRATOS_ERROR` |
| `Info()`, `PrintInfo()`, `PrintData()`, `save` / `load` | all | Standard Kratos boilerplate (the laws are serializable) |
| `static constexpr double ZeroTolerance` | `FrictionalLaw` | `std::numeric_limits<double>::epsilon()` |

`FrictionalLawWithDerivative` ([header](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_frictional_laws/frictional_law_with_derivative.h)) introduces the template parameters of the mortar conditions and two typedefs that tie it to the derivative machinery of the Kratos core (`kratos/includes/mortar_classes.h`): `DerivativeDataType = DerivativeDataFrictional<TDim, TNumNodes, TNumNodesMaster>` (current and reference coordinates `u1`, `X1`, `u2`, `X2`, normal derivatives `DeltaNormalSlave`, previous multipliers) and `MortarConditionMatrices = MortarOperatorWithDerivatives<TDim, TNumNodes, true, TNumNodesMaster>` (the mortar operators $$\mathbf{D}$$, $$\mathbf{M}$$ and their directional derivatives `DeltaDOperator[i]`, `DeltaMOperator[i]`). See [Linearisation and derivatives](../Theory/Linearisation_And_Derivatives.html) for what these objects contain.

### Coulomb

[`coulomb_frictional_law.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_frictional_laws/coulomb_frictional_law.cpp). Coulomb's law (thesis eq. 4.45a, $$\phi_{co} := \Vert\mathbf{t}_\tau\Vert - \mu\Vert p_n\Vert \le 0$$) bounds the tangential traction by the friction coefficient times the normal pressure. In the augmented Lagrangian formulation the normal pressure of a slave node is the *augmented* pressure $$\bar{\lambda}_n = k\,\mathbf{n}\cdot\boldsymbol{\lambda} + \varepsilon\tilde{g}_n$$ (negative in compression), so the threshold implemented by `GetThresholdValue` is

<p align="center">$$\mathscr{F}_{ALM} = -\mu\,\bar{\lambda}_n = -\mu\,\big(k\,\mathbf{n}\cdot\boldsymbol{\lambda} + \varepsilon\,\tilde{g}_n\big)$$</p>

(thesis eq. 4.76). The code reads $$\bar{\lambda}_n$$ from the nodal non-historical variable `AUGMENTED_NORMAL_CONTACT_PRESSURE`, which `ActiveSetUtilities` stores at every active-set check (see [Strategies and convergence criteria](Strategies_And_Convergence_Criteria.html)); a block that recomputes it by hand from `SCALE_FACTOR`, `INITIAL_PENALTY`, `VECTOR_LAGRANGE_MULTIPLIER`, `NORMAL` and `WEIGHTED_GAP` is kept in the source as a comment (`coulomb_frictional_law.cpp:29-40`).

`GetDerivativeThresholdValue` linearizes $$\mathscr{F}_{ALM}$$ with respect to the DoFs of the pair. Let $$n_u = d\,(n_s + n_m)$$ be the number of displacement DoFs of the pair (slave nodes first, master nodes after). For `IndexDerivative` $$\ge n_u$$ the derivative is with respect to a multiplier component; it is non-zero only for the multiplier of the node itself (`aux_index / TDim == IndexNode`) and equals $$-\mu\,k\,n_j$$ with $$j$$ the component (`aux_index % TDim`). For `IndexDerivative` $$\lt n_u$$ the derivative is with respect to a displacement and only the gap term contributes:

<p align="center">$$\frac{\partial\mathscr{F}_{ALM}}{\partial u} = -\mu\,\varepsilon\,\frac{\partial\tilde{g}_n}{\partial u}, \qquad \frac{\partial\tilde{g}_n}{\partial u} = -\mathbf{n}\cdot\Big[\big(\mathbf{M}\,\delta\mathbf{x}_1 - \mathbf{D}\,\delta\mathbf{x}_2\big) + \big(\Delta\mathbf{M}\,\mathbf{x}_1 - \Delta\mathbf{D}\,\mathbf{x}_2\big)\Big]_{\text{row } IndexNode} - \Delta\mathbf{n}\cdot\big[\mathbf{M}\,\mathbf{x}_1 - \mathbf{D}\,\mathbf{x}_2\big]_{\text{row } IndexNode}$$</p>

where $$\delta\mathbf{x}_1$$ / $$\delta\mathbf{x}_2$$ are unit perturbations of the slave / master coordinate selected by `IndexDerivative` (`Deltax1`, `Deltax2`), $$\mathbf{x}_1 = \mathbf{u}_1 + \mathbf{X}_1$$, $$\mathbf{x}_2 = \mathbf{u}_2 + \mathbf{X}_2$$ are the current coordinates from `rDerivativeData`, $$\Delta\mathbf{D}$$, $$\Delta\mathbf{M}$$ are `DeltaDOperator[IndexDerivative]`, `DeltaMOperator[IndexDerivative]`, and the last term (normal variation, `DeltaNormalSlave[IndexDerivative]`) is only added when `TNormalVariation` is `true` and the derivative is taken with respect to a slave node. Note: the source uses the mortar operators with swapped roles in the two products (`prod(MOperator, Deltax1) - prod(DOperator, Deltax2)` and `prod(DeltaMOperator, x1) - prod(DeltaDOperator, x2)`), while the weighted gap of the conditions is $$\tilde{g}_n = -\mathbf{n}\cdot(\mathbf{D}\mathbf{x}_1 - \mathbf{M}\mathbf{x}_2)$$; since the class is not called by the conditions at present (see status below) this has no effect on the results.

### Tresca

[`tresca_frictional_law.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_frictional_laws/tresca_frictional_law.cpp). Tresca's model uses a threshold independent of the normal pressure ("simpler because does not depend on contact normal pressure, it is just a threshold value", thesis footnote 8 on p. 113):

<p align="center">$$\mathscr{F}_{Tresca} = \tau_{max}$$</p>

`GetThresholdValue` looks up the application variable `TRESCA_FRICTION_THRESHOLD` (a `double`) with the same precedence as the friction coefficient (condition properties, node, `ProcessInfo`) and returns 0.0 when it is not found; `GetDerivativeThresholdValue` returns 0.0 because the threshold is constant.

### Template instantiations and Python names

Each templated law is explicitly instantiated at the bottom of its `.cpp` for the ten combinations used by the conditions, `<2,2,·,2>`, `<3,3,·,3>`, `<3,4,·,4>`, `<3,3,·,4>`, `<3,4,·,3>` with `TNormalVariation` `false` and `true`. [`add_custom_frictional_laws_to_python.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_python/add_custom_frictional_laws_to_python.cpp) exposes the base `FrictionalLaw` (default constructor, `GetFrictionCoefficient`, `GetThresholdValue`) and, through `RegisterFrictionalLaws<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>(m, rEndName)`, the classes `FrictionalLaw<suffix>`, `TrescaFrictionalLaw<suffix>` and `CoulombFrictionalLaw<suffix>` (default constructors only) for:

| Suffix | Template arguments | Geometry pair |
|---|---|---|
| `2D2N` / `2D2NNV` | `<2, 2, false/true, 2>` | line – line |
| `3D3N` / `3D3NNV` | `<3, 3, false/true, 3>` | triangle – triangle |
| `3D4N` / `3D4NNV` | `<3, 4, false/true, 4>` | quadrilateral – quadrilateral |
| `3D3N4N` / `3D3N4NNV` | `<3, 3, false/true, 4>` | triangle slave – quadrilateral master |
| `3D4N3N` / `3D4N3NNV` | `<3, 4, false/true, 3>` | quadrilateral slave – triangle master |

`NV` stands for *normal variation* (`TNormalVariation = true`), exactly as in the condition names `ALMNVFrictionalMortarContactCondition…` (see [Conditions](Conditions.html)). Example: `KratosMultiphysics.ContactStructuralMechanicsApplication.CoulombFrictionalLaw3D3NNV()`.

### Status of the integration

The frictional-law layer is complete as a set of classes but is **not yet the source of the friction coefficient used by the conditions**. The top-level `README.md` lists "Frictional laws (WIP) in order to consider different types of frictional behaviour". Concretely:

- The application variable `FRICTIONAL_LAW` (type `FrictionalLaw::Pointer`) is declared, created and registered (`contact_structural_mechanics_application_variables.{h,cpp}`, `contact_structural_mechanics_application.cpp:178`), but it is not exposed to Python and **no C++ or Python code reads or writes it**.
- The frictional conditions obtain the friction coefficient nodally: `GetFrictionCoefficient()` in `ALM_frictional_mortar_contact_condition.h:456-467` and `penalty_frictional_mortar_contact_condition.h` fills a vector with `r_geometry[i_node].GetValue(FRICTION_COEFFICIENT)` and carries the comment `// TODO: Define the "CL" or friction law to compute this`. The nodal values are set by `ALMFastInit` (`alm_fast_init_process.cpp:76-112`): for a frictional problem it zeroes `FRICTION_COEFFICIENT` and `NODAL_AREA` on the interface nodes, adds the `FRICTION_COEFFICIENT` of the properties of every condition to its nodes (warning "Friction coefficient not defined, zero will be considered" when a property lacks it) together with a counter in `NODAL_AREA`, and divides, so a node shared by pairs with different coefficients gets their arithmetic mean. The properties themselves receive the coefficient from the `friction_coefficients` dictionary of the contact process (`alm_contact_process.py:463-476`, `mpc_contact_process.py:399-402`), with a warning if the property already had one.
- The active-set utilities (`ActiveSetUtilities::ComputeALMFrictionalActiveSet`, `ComputePenaltyFrictionalActiveSet`) and `MPCContactCriteria` also read the nodal `FRICTION_COEFFICIENT` and evaluate the Coulomb threshold $$-\mu\bar{\lambda}_n$$ inline, without going through a `FrictionalLaw` object.
- The contact processes accept the key `"frictional_law" : "Coulomb"` (`alm_contact_process.py:70`, `penalty_contact_process.py:52`, `explicit_penalty_contact_process.py:52`, `mpc_contact_process.py:69`) and store it in `self.frictional_law`, but the value is never used afterwards; only the Coulomb behavior hard-coded in the generated conditions and in the utilities is available. Note: setting `"frictional_law" : "Tresca"` has no effect.
- `TRESCA_FRICTION_THRESHOLD` is registered and exposed to Python, so it can be assigned to properties or nodes, but it is only read by `TrescaFrictionalLaw::GetThresholdValue`.

The generated frictional conditions therefore implement Coulomb friction with a nodal, possibly heterogeneous, coefficient; the `FrictionalLaw` classes are the intended extension point for the day the conditions delegate the threshold and its derivative (which is what the `GetDerivativeThresholdValue` signature, with the `DerivativeData` and `MortarOperatorWithDerivatives` arguments, is prepared for). The friction cone and the regularized Coulomb law are discussed in [Frictional contact](../Theory/Frictional_Contact.html) (thesis §4.2.5 and Fig. 4.10).

## The MPC contact route

<p align="center"><img src="images/csma_mpc_contact_flow.svg" alt="The MPC contact route: mpc_contact_process.py, MPCContactSearchProcess, MPCMortarContactCondition with ContactMasterSlaveConstraint, MPC strategy and MPCContactCriteria" width="1000"/></p>
<p align="center"><em>Figure: the MPC contact route from the Python settings to the eliminated system (thesis App. D.5).</em></p>

In the multiplier and penalty formulations the contact condition assembles a local system. In the MPC route the same mortar operators are used to write the **slave displacement as a linear function of the master displacements**,

<p align="center">$$\mathbf{u}_{s} = \mathbf{D}^{-1}\mathbf{M}\,\mathbf{u}_{m} + \mathbf{g}, \qquad \mathbf{D}^{-1} \text{ trivial for dual multipliers},$$</p>

and to hand this relation to the Kratos `MasterSlaveConstraint` machinery, which eliminates the slave DoFs at the builder-and-solver level (thesis App. D.5, "MultiPoint Constraint (Master-Slave elimination method)", whose §D.5.3 discusses its applicability to contact). The contact inequality is then handled by activating and deactivating constraints with an active-set check based on the reactions. The route has no multipliers, no penalty and no condensation, at the price of a coarser treatment of the non-linearity (the relation matrix is frozen during the Newton iterations unless `update_each_nl_iteration` is set) and of a heuristic reaction-based active set. The classes involved are:

| Layer | Class | Where documented |
|---|---|---|
| Python process | `MPCContactProcess` (`mpc_contact_process.py`), derived from `SearchBaseProcess` | [Processes](Processes.html), [Contact process settings reference](../Usage/Contact_Process_Settings_Reference.html) |
| Python solver | `MPCContactStaticSolver`, `MPCContactImplicitMechanicalSolver` (`mpc_contact_structural_mechanics_{static,implicit_dynamic}_solver.py`), selected when `solver_settings` has `mpc_contact_settings` | [Solver settings reference](../Usage/Solver_Settings_Reference.html) |
| Search | `MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>` and `MPCContactSearchWrapperProcess` | [Search pipeline and bounding volumes](../Contact_Search/Search_Pipeline_And_Bounding_Volumes.html) |
| Condition | `MPCMortarContactCondition<TDim, TNumNodes, TNumNodesMaster>` (`MPCMortarContactCondition2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N`) | [Conditions](Conditions.html), this page |
| Constraint | `ContactMasterSlaveConstraint` | this page |
| Strategy and criterion | `ResidualBasedNewtonRaphsonMPCContactStrategy`, `MPCContactCriteria` | [Strategies and convergence criteria](Strategies_And_Convergence_Criteria.html) |
| Builder-and-solver | Core `ResidualBasedBlockBuilderAndSolver` (default) or `ResidualBasedEliminationBuilderAndSolverWithConstraints`; the contact-specific `ContactResidualBasedEliminationBuilderAndSolverWithConstraints` is selected by the mortar-based solvers | [Builder and solvers and linear solvers](Builder_And_Solvers_And_Linear_Solvers.html) |

### `ContactMasterSlaveConstraint`

[`contact_master_slave_constraint.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_master_slave_constraints/contact_master_slave_constraint.h) / [`.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_master_slave_constraints/contact_master_slave_constraint.cpp). The only class of `custom_master_slave_constraints/`; it derives from the core `LinearMasterSlaveConstraint` (`BaseType`, whose own base is `MasterSlaveConstraint`) and is registered as `"ContactMasterSlaveConstraint"` in `KratosContactStructuralMechanicsApplication::Register()` (`KRATOS_REGISTER_CONSTRAINT`). A linear master–slave constraint stores a list of slave DoFs, a list of master DoFs, a relation matrix $$\mathbf{T}$$ and a constant vector $$\mathbf{c}$$, meaning $$\mathbf{u}_{slave} = \mathbf{T}\,\mathbf{u}_{master} + \mathbf{c}$$. The contact version adds nothing to the data; its purpose is to be a *distinct type* that the search can create, the condition can update and the criterion can (de)activate.

| Member | Behavior |
|---|---|
| `ContactMasterSlaveConstraint(IndexType Id = 0)` | Empty constraint, the one created by the search |
| `ContactMasterSlaveConstraint(Id, rMasterDofsVector, rSlaveDofsVector, rRelationMatrix, rConstantVector)` | Full constructor forwarding to `LinearMasterSlaveConstraint` |
| `ContactMasterSlaveConstraint(Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant)` | Single-component constructor; **always throws** `KRATOS_ERROR` ("Please don't use this constructor. A components variable is expected") |
| `Create(...)` (two overloads) | Return `Kratos::make_shared<ContactMasterSlaveConstraint>` with the same arguments |
| `FinalizeNonLinearIteration(const ProcessInfo&) override` | **Empty body** (`contact_master_slave_constraint.cpp:115-118`); the update of the relation is driven by the condition, not by the constraint |
| `GetInfo()` | Returns `"This is contact MPC !"`; `PrintInfo` writes the class name and id |
| `save` / `load` | Serialize through the base class |

The relation of a constraint is (re)written from outside with the inherited `SetDofList(rSlaveDofsVector, rMasterDofsVector, rCurrentProcessInfo)` and `SetLocalSystem(rRelationMatrix, rConstantVector, rCurrentProcessInfo)`; its participation in the system is controlled with the `ACTIVE` flag, which the builder-and-solver honors.

### Creation and attachment: `MPCContactSearchProcess`

`MPCContactSearchProcess` ([`mpc_contact_search_process.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_processes/mpc_contact_search_process.cpp)) derives from `BaseContactSearchProcess` and runs the same KD-tree / OBB search as the other formulations (see [Search pipeline and bounding volumes](../Contact_Search/Search_Pipeline_And_Bounding_Volumes.html)); its overrides deal with the constraints:

- `AddPairing` (lines 113-145) creates the pair condition through the base class (the condition name is `MPCMortarContact` + `Condition` + dimension/node suffix, from `MPCContactProcess._get_condition_name`), then a `ContactMasterSlaveConstraint` with id `GetMaximumConstraintsIds() + 1`, sets it `ACTIVE`, initializes it, adds it to the computing model part and **stores its pointer in the condition** with `p_cond->SetValue(CONSTRAINT_POINTER, p_new_const)`. The condition inherits the problem type from the main model part flags: `SLIP` if the main model part `Is(SLIP)` (frictional), `RIGID` if it `Is(RIGID)` (mesh tying through constraints), nothing for frictionless contact. Finally `p_cond->Initialize(r_process_info)` builds the first relation (see below).
- `CheckContactModelParts` (lines 48-108) clones the constraints of the `Contact` sub-model-part that carry the `MARKER` flag with new ids (so that a constraint that survives a remeshing is not shared between the old and new sub-model-parts) and marks the others, removing the `TO_ERASE` ones.
- `CleanModelPart` (lines 149-160) flags all constraints of the given model part `TO_ERASE` and removes them from all levels; `ResetContactOperators` (lines 184-260) does the same for the constraints of the `ComputingContact(Sub<id>)` model part after a remeshing (`MODIFIED`), or removes only the constraints of the pairs that disappeared, using the `INDEX_MAP` of every master node.

`CONSTRAINT_POINTER` is an application variable of type `MasterSlaveConstraint::Pointer` registered in `contact_structural_mechanics_application_variables.h`; it is not exposed to Python (see the [Variables and flags reference](Variables_And_Flags_Reference.html)).

### Update: `MPCMortarContactCondition`

[`mpc_mortar_contact_condition.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_conditions/mpc_mortar_contact_condition.cpp). The condition derives directly from `PairedCondition` (not from `MortarContactCondition`), has `MatrixSize = TDim * (TNumNodes + TNumNodesMaster)` and **assembles nothing**: `CalculateLeftHandSide` and `CalculateRightHandSide` (lines 351-410) resize and zero the local system, `CalculateMassMatrix` / `CalculateDampingMatrix` likewise. Its role is to keep the constraint attached through `CONSTRAINT_POINTER` up to date with the current geometry. The same update block appears three times:

| Hook | When the constraint is rebuilt |
|---|---|
| `Initialize` (lines 78-121) | Always, right after the search created the pair (with an empty auxiliary `ProcessInfo`) |
| `InitializeSolutionStep` (lines 125-172) | At the beginning of every time step unless the condition `Is(BLOCKED)`. `MPCContactProcess.ExecuteFinalizeSolutionStep` sets `BLOCKED` on all conditions when `update_condition_relation_step` is `false` (the default), so by default the relation is computed once, at the step in which the pair is created; with `update_condition_relation_step = true` it follows the deformation step by step. For frictional conditions the previous mortar operators are stored here the first time (`mPreviousMortarOperatorsInitialized`) |
| `InitializeNonLinearIteration` (lines 176-216) | At every Newton iteration, only if the condition `Is(INTERACTION)`, which `ResidualBasedNewtonRaphsonMPCContactStrategy::AuxiliarySolveSolutionStep` sets from `update_each_nl_iteration` |

The block itself:

1. `MortarExplicitContributionUtilities<TDim, TNumNodes, FRICTIONLESS_PENALTY, false, TNumNodesMaster>::ComputePreviousMortarOperators(this, rCurrentProcessInfo, mortar_operators, integration_order, false)` integrates the mortar operators $$\mathbf{D}$$ and $$\mathbf{M}$$ of the pair on the current configuration (`INTEGRATION_ORDER_CONTACT` from the properties, default 2) and returns whether **dual** shape functions were used (`dual_LM`).
2. A relation matrix of size `TDim * TNumNodes × TDim * TNumNodesMaster` and a constant vector of size `TDim * TNumNodes` are zeroed and filled by `UpdateConstraintFrictional` (condition `Is(SLIP)`), `UpdateConstraintTying` (condition `Is(RIGID)`) or `UpdateConstraintFrictionless` (otherwise).
3. `ConstraintDofDatabaseUpdate` prunes the relation and rebuilds the DoF lists of the constraint.
4. `p_const->SetLocalSystem(relation_matrix, constant_vector, rCurrentProcessInfo)` stores the new relation.

**`UpdateConstraintFrictionless`** (lines 593-710). With $$\mathbf{D}$$ and $$\mathbf{M}$$ from the mortar operators, $$\mathbf{D}^{-1}$$ is the trivial diagonal inverse when `dual_LM` is `true` and a full `MathUtils::InvertMatrix` otherwise; the product $$\mathbf{D}^{-1}\mathbf{M}$$ (`D_inv_M`, `TNumNodes × TNumNodesMaster`) is the mortar projection of master onto slave nodes. Any row of $$\mathbf{D}^{-1}\mathbf{M}$$ with an entry below $$-10^{-8}$$ is zeroed entirely (a negative weight means the slave node projects outside the master, "not in the same section"). Then, for every **active** slave node $$i$$ (`IsActive()`), with $$w_i = 1/$$`NODAL_PAUX` (the inverse of the number of pair conditions sharing the node, computed by the strategy's `ComputeNodalWeights`, so that a node shared by several pairs receives the average of their relations) and $$\mathbf{n}_i$$ the nodal normal, the *normal* components are coupled:

<p align="center">$$T_{(i,a),(j,b)} = w_i\, n_{i,a}\, n_{i,b}\,\big(\mathbf{D}^{-1}\mathbf{M}\big)_{ij}, \qquad a, b = 1..d,\ j = 1..n_m$$</p>

(entries with $$\vert T\vert \le 10^{-12}$$ are not written). This is a *normal-only* constraint: the slave node is free to move tangentially, which is what frictionless contact requires. The constant vector carries the gap that must be closed. The nodal gap is $$g_i = \tilde{g}_{n,i}/$$`NODAL_AREA` (the weighted gap normalized by the nodal area, both accumulated by `ComputeExplicitContributionConditions`), corrected by the *previous* displacements so that the constraint is written in terms of total displacements: $$g_i \leftarrow g_i + \mathbf{n}_i\cdot\big[\mathbf{u}_1^{n-1} - \mathbf{D}^{-1}\mathbf{M}\,\mathbf{u}_2^{n-1}\big]_i$$ (`u1_0`, `u2_0` are `DISPLACEMENT` at buffer position 1), and

<p align="center">$$c_{(i,a)} = w_i\, g_i\, n_{i,a}.$$</p>

A variant that computes the gap from the positions instead of the previous displacements is kept as a comment. Inactive slave nodes get zero rows, which `ConstraintDofDatabaseUpdate` then removes.

**`UpdateConstraintFrictional`** (lines 716-826). Same $$\mathbf{D}^{-1}\mathbf{M}$$ and same constant vector; the relation matrix depends on the frictional state of each active slave node: a `SLIP` node receives the normal-only coupling above (comment: `// SLIP // TODO: Add nodal forces`), a stick node receives the full vector coupling

<p align="center">$$T_{(i,a),(j,a)} = w_i\,\big(\mathbf{D}^{-1}\mathbf{M}\big)_{ij}, \qquad a = 1..d$$</p>

(comment: `// STICK // TODO: ADD the contribution of slip to constant vector`), which glues the slave node to the master surface in all directions. Note: the frictional MPC route therefore treats stick exactly and slip as frictionless (no tangential force is imposed on slipping nodes); the `SLIP` flag is toggled by `MPCContactCriteria` from the mapped reactions and the nodal `FRICTION_COEFFICIENT`, and the `FrictionalLaw` classes are not involved.

**`UpdateConstraintTying`** (lines 832-925). Used when the main model part `Is(RIGID)`, that is when `mpc_contact_settings.contact_type` / the process `contact_type` contains `MeshTying`: the full vector coupling of the stick case is written for every active slave node without looking at `SLIP`, and the same gap-based constant vector is built. Since the active set of `MPCContactCriteria` still applies (`RIGID` is only checked when the problem is not frictional), the tying constraint is released for nodes whose reaction is in traction; this is the "tying with tension check" of the flowchart above.

**`ConstraintDofDatabaseUpdate`** (lines 466-587). The relation matrix produced above has many zero rows (inactive slave nodes, tangential components in the frictionless case) and possibly zero columns (master nodes with no projection). A constraint with zero rows would declare slave DoFs that are not actually constrained, so the method computes the sum of absolute values of every row and column, keeps only the rows (`slave_dofs_OK`) and columns (`master_dofs_OK`) whose sum exceeds $$10^{-4}$$, rebuilds a dense relation matrix and constant vector of the reduced size, and builds the matching slave and master `Dof` lists (`DISPLACEMENT_X`, `_Y` and, in 3D, `_Z` of the parent and paired geometries, in the same order as the rows/columns). Finally `p_const->SetDofList(slave_dof_vector, master_dof_vector, rCurrentProcessInfo)`. When nothing was pruned the full lists are used. A commented-out block shows an earlier design in which the slave DoFs were also listed as masters (identity rows for inactive nodes), which was abandoned in favor of pruning.

Two further details of the condition: `AddExplicitContribution` (lines 447-462) calls `MortarExplicitContributionUtilities::AddExplicitContributionOfMortarCondition` (or the frictional variant with the previous operators) so that `WEIGHTED_GAP`, `WEIGHTED_SLIP` and `NODAL_AREA` are accumulated exactly as in the other formulations — this is what `MPCContactCriteria::PreCriteria` / `PostCriteria` and the strategy's `Predict` rely on; and `ComputePreviousMortarOperators` (lines 935-940) stores $$\mathbf{D}^{n-1}$$, $$\mathbf{M}^{n-1}$$ for the slip computation of frictional problems (called in `InitializeSolutionStep` and `FinalizeSolutionStep`).

### Activation, deactivation and the tension check

`MPCContactCriteria::PostCriteria` (see [Strategies and convergence criteria](Strategies_And_Convergence_Criteria.html#mpccontactcriteria)) maps the master `REACTION` onto the slave nodes, computes for every slave node the normal contact pressure $$p_n = (\mathbf{R}\cdot\mathbf{n})/$$`NODAL_MAUX` and the normalized gap $$g = \tilde{g}_n/$$`NODAL_AREA`, and declares the node **active** when

<p align="center">$$p_n \lt -\,\texttt{REACTION\_CHECK\_STIFFNESS\_FACTOR}\cdot E \quad\text{or}\quad g \lt 0$$</p>

where $$E$$ is the `YOUNG_MODULUS` of the first element properties. The factor (default $$10^{-10}$$ from `mpc_contact_settings`, written to `ProcessInfo[REACTION_CHECK_STIFFNESS_FACTOR]` by `MPCContactProcess._initialize_process_info`) turns the reaction check into a relative one: a slave node stays active while its reaction pushes into the master with a pressure larger than a tiny fraction of the stiffness, or while it penetrates; it is released when the constraint would have to *pull* the node (tension). A pair condition whose slave nodes are all inactive is set inactive together with its constraint (`p_const->Set(ACTIVE, false)`), so that the builder-and-solver ignores it in the next assembly; the strategy re-runs `ComputeNodalWeights` before every iteration and, with `update_each_nl_iteration`, `SetUpDofSet` / `SetUpSystem`, because the number of eliminated DoFs changes with the active set. Nodes are re-activated by the same test in the next `PostCriteria`; the loop converges when no node changes status (and, for frictional problems, when no node changes between stick and slip).

### End to end

1. **Settings.** `solver_settings.mpc_contact_settings` (`contact_type`, `simplified_semi_smooth_newton`, `inner_loop_iterations` = 10, `update_each_nl_iteration`, `enforce_ntn`) selects `MPCContactStaticSolver` / `MPCContactImplicitMechanicalSolver` through the structural solver wrapper; `processes.contact_process_list` holds an `mpc_contact_process` (`MPCContactProcess`) with the same pairing dictionaries as the ALM process plus `reaction_check_stiffness_factor`, `tangent_factor`, `zero_tolerance_factor`, `frictional_law` (accepted, unused) and `update_condition_relation_step`.
2. **Set-up (Python).** `MPCContactProcess.ExecuteInitialize` creates the `Contact`, `ContactSub<key>`, `MasterSubModelPart<key>` and `SlaveSubModelPart<key>` model parts, sets `SLIP` (frictional) or `RIGID` (mesh tying) on the main model part, writes `REACTION_CHECK_STIFFNESS_FACTOR`, `TANGENT_FACTOR`, `ACTIVE_SET_CONVERGED` and the search parameters to the `ProcessInfo`, runs `ALMFastInit` (normals, `WEIGHTED_GAP`, friction coefficients) and creates one `MPCContactSearchProcess` per pair (`_create_main_search`).
3. **Search (C++).** At every step (`database_step_update`) `MPCContactSearchProcess` clears the stale pairs and constraints and, for every new slave–master pair, creates an `MPCMortarContactCondition` in `ComputingContact` and a `ContactMasterSlaveConstraint` linked through `CONSTRAINT_POINTER`; `Initialize` of the condition writes the first relation.
4. **Solve.** `ResidualBasedNewtonRaphsonMPCContactStrategy` predicts, computes the nodal weights, and runs the Newton loop; `MPCMortarContactCondition::InitializeSolutionStep` / `InitializeNonLinearIteration` refresh the relations, the builder-and-solver eliminates the slave DoFs, and `MPCContactCriteria` (both the strategy's internal one and the one in the user `AndCriteria`) updates the active set from the reactions until it is stable.
5. **Finalize.** `MPCContactProcess.ExecuteFinalizeSolutionStep` blocks the conditions (unless `update_condition_relation_step`), resets the `SLIP` flags of inactive nodes for post-processing and the search removes the pairs that are no longer close.

The tests of this route are in `tests/mpc_contact_tests/` (`beam_contact_static_test`, `beam_contact_static_with_friction_test`, `plate_test` with `contact_type` `MeshTying`, `3D_multi_contact_test`, …, see the [Test suite reference](../Validation/Test_Suite_Reference.html)); they use the default `block` builder-and-solver of the structural solver. For the theory of the master–slave elimination and its comparison with penalty and Lagrange multipliers see [Constrained optimisation methods](../Theory/Constrained_Optimisation_Methods.html); for the mortar operators $$\mathbf{D}$$ and $$\mathbf{M}$$ that make $$\mathbf{D}^{-1}\mathbf{M}$$ cheap see [Mortar integration and dual Lagrange multipliers](../Theory/Mortar_Integration_And_Dual_Lagrange_Multipliers.html).

## Quick reference

| Item | Value |
|---|---|
| Frictional law classes | `FrictionalLaw`, `FrictionalLawWithDerivative<…>`, `CoulombFrictionalLaw<…>`, `TrescaFrictionalLaw<…>` (10 instantiations each) |
| Python names | `FrictionalLaw`, `FrictionalLaw<S>`, `CoulombFrictionalLaw<S>`, `TrescaFrictionalLaw<S>` with `S` in `2D2N, 3D3N, 3D4N, 3D3N4N, 3D4N3N` and the `NV` versions |
| Coulomb threshold | $$-\mu\,\bar{\lambda}_n$$, $$\bar{\lambda}_n$$ read from `AUGMENTED_NORMAL_CONTACT_PRESSURE` |
| Tresca threshold | `TRESCA_FRICTION_THRESHOLD` (properties → node → `ProcessInfo`) |
| Friction coefficient lookup | `FRICTION_COEFFICIENT` (properties → node → `ProcessInfo`) in the law; nodal value averaged by `ALMFastInit` in the conditions |
| Unused hooks | `FRICTIONAL_LAW` variable, `frictional_law` JSON key |
| Constraint class | `ContactMasterSlaveConstraint` (registered `"ContactMasterSlaveConstraint"`), empty `FinalizeNonLinearIteration` |
| Attachment | `condition.SetValue(CONSTRAINT_POINTER, constraint)` in `MPCContactSearchProcess::AddPairing` |
| Relation | $$\mathbf{T} = w\,(\mathbf{n}\otimes\mathbf{n})\,\mathbf{D}^{-1}\mathbf{M}$$ (frictionless, slip) or $$\mathbf{T} = w\,\mathbf{D}^{-1}\mathbf{M}\otimes\mathbf{I}$$ (stick, tying); $$\mathbf{c} = w\,g\,\mathbf{n}$$ |
| Pruning thresholds | Row/column sum $$\gt 10^{-4}$$ kept; entries $$\gt 10^{-12}$$ written; rows with an entry $$\lt -10^{-8}$$ in $$\mathbf{D}^{-1}\mathbf{M}$$ zeroed |
| Active-set test | $$p_n \lt -\texttt{REACTION\_CHECK\_STIFFNESS\_FACTOR}\cdot E$$ or $$g \lt 0$$ |
| Relation update frequency | Once per pair (default), per step (`update_condition_relation_step`), per iteration (`update_each_nl_iteration`) |
