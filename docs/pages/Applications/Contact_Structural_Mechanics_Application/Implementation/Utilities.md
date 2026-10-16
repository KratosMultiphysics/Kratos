---
title: Utilities
keywords: contact, utilities, ContactUtilities, ActiveSetUtilities, SelfContactUtilities, DerivativesUtilities, MortarExplicitContributionUtilities, InterfacePreprocessCondition, ProcessFactoryUtility, MortarUtilities, ExactMortarIntegrationUtility, mortar_classes
tags: [contact, implementation, utilities, active set, mortar, derivatives]
sidebar: contact_structural_mechanics_application
summary: Reference of the utility classes and namespaces of the ContactStructuralMechanicsApplication (ContactUtilities, ActiveSetUtilities, SelfContactUtilities, DerivativesUtilities, MortarExplicitContributionUtilities, InterfacePreprocessCondition, ProcessFactoryUtility, logging macros), their public API and Python exposure, the Kratos-core mortar utilities the application depends on, and the mapping between utilities and their C++ unit tests.
---

> **Sources.** Thesis §4.3.3.5 (active set strategy, pp. 105–108), §4.4.4 (gap and pairing, pp. 130–137), §4.6 (directional derivatives, pp. 141–160); code: `custom_utilities/contact_utilities.{h,cpp}`, `custom_utilities/active_set_utilities.{h,cpp}`, `custom_utilities/self_contact_utilities.{h,cpp}`, `custom_utilities/derivatives_utilities.{h,cpp}`, `custom_utilities/mortar_explicit_contribution_utilities.{h,cpp}`, `custom_utilities/interface_preprocess.{h,cpp}`, `custom_utilities/logging_settings.hpp`, `custom_python/process_factory_utility.{h,cpp}`, `custom_python/add_custom_utilities_to_python.cpp`; Kratos core `kratos/utilities/mortar_utilities.h`, `kratos/utilities/exact_mortar_segmentation_utility.h`, `kratos/includes/mortar_classes.h`; tests `tests/cpp_tests/utilities/*.cpp`.

The utilities are the stateless helpers shared by the [conditions](Conditions.html), the [processes](Processes.html) and the [strategies and convergence criteria](Strategies_And_Convergence_Criteria.html). Unlike processes they have no life cycle: they are static classes or namespaces of free functions that operate on a model part, a geometry or a condition. This page lists their public API, states which functions are reachable from Python and points to the unit tests that pin down their behavior. The mathematical background of the active set is in [Frictionless contact](../Theory/Frictionless_Contact.html) and [Frictional contact](../Theory/Frictional_Contact.html), the derivatives in [Linearisation and derivatives](../Theory/Linearisation_And_Derivatives.html) and the mortar operators in [Mortar integration and dual Lagrange multipliers](../Theory/Mortar_Integration_And_Dual_Lagrange_Multipliers.html).

## Overview

| Component | File(s) | Kind | Python exposure | Main clients |
|---|---|---|---|---|
| `ContactUtilities` | `custom_utilities/contact_utilities.{h,cpp}` | static class | `CSMA.ContactUtilities` (8 of 14 functions) | search processes, criteria, strategies, Python processes |
| `ActiveSetUtilities` | `custom_utilities/active_set_utilities.{h,cpp}` | namespace | submodule `CSMA.ActiveSetUtilities` (2 of 5 functions) | mortar convergence criteria, `ExplicitPenaltyContactProcess` |
| `SelfContactUtilities` | `custom_utilities/self_contact_utilities.{h,cpp}` | namespace | submodule `CSMA.SelfContactUtilities` (3 of 3) | `BaseContactSearchProcess` |
| `DerivativesUtilities<TDim,TNumNodes,TFrictional,TNormalVariation,TNumNodesMaster>` | `custom_utilities/derivatives_utilities.{h,cpp}` | static template class | none | mortar contact conditions |
| `MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional,TNormalVariation,TNumNodesMaster>`, `AuxiliaryOperationsUtilities` | `custom_utilities/mortar_explicit_contribution_utilities.{h,cpp}` | static template class + namespace | none | `MortarContactCondition::AddExplicitContribution`, penalty conditions |
| `InterfacePreprocessCondition` | `custom_utilities/interface_preprocess.{h,cpp}` | class | `CSMA.InterfacePreprocessCondition` | `SearchBaseProcess._interface_preprocess` |
| `ProcessFactoryUtility` | `custom_python/process_factory_utility.{h,cpp}` | class | `CSMA.ProcessFactoryUtility` | contact solvers → contact strategies |
| logging macros | `custom_utilities/logging_settings.hpp` | preprocessor macros | none | debugging of conditions and the mixed linear solver |

`CSMA` stands for `KratosMultiphysics.ContactStructuralMechanicsApplication`. The bindings live in [`add_custom_utilities_to_python.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_python/add_custom_utilities_to_python.cpp).

## `ContactUtilities`

[`contact_utilities.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_utilities/contact_utilities.h) groups generic helpers on model parts and geometries. All functions are `static`; from Python they are called on the class (`CSMA.ContactUtilities.CheckActivity(model_part)`), a default constructor being exposed only for convenience.

| Function | Signature | Python | Purpose | Used by |
|---|---|---|---|---|
| `CalculateRelativeSizeMesh` | `double (ModelPart&)` | yes | `CalculateMaxNodalH / CalculateMinimalNodalH` | `SearchBaseProcess` (`adapt_search`) |
| `CalculateMaxNodalH` | `double (ModelPart&)` | yes | Maximum historical `NODAL_H` (parallel `MaxReduction`) | search (`bounding_box_factor` scaling) |
| `CalculateMeanNodalH` | `double (ModelPart&)` | yes | Mean `NODAL_H` | `AdvancedContactSearchProcess` (`DISTANCE_THRESHOLD`), `ExplicitPenaltyContactProcess` (`MAX_GAP_THRESHOLD`) |
| `CalculateMinimalNodalH` | `double (ModelPart&)` | yes | Minimum `NODAL_H` | `CalculateRelativeSizeMesh` |
| `ScaleNode<TPointType>` | `void (TPointType&, const array_1d<double,3>& rNormal, const double LengthSearch)` | no | Moves a point along its normal by `LengthSearch` (enlarges the search box) | `BaseContactSearchProcess` (`InBox` search) |
| `DistancePoints` | `double (const CoordinatesArrayType&, const CoordinatesArrayType&)` | no | Euclidean distance | search |
| `ComputeStepJump` | `void (ModelPart&, const double DeltaTime, const bool HalfJump = true)` | no | Writes in `DELTA_COORDINATES` the predicted displacement $$\Delta t\,\mathbf{v} (+ \tfrac{1}{2}\Delta t^2 \mathbf{a})$$ of every node (half step when `HalfJump`) | `BaseContactSearchProcess::UpdatePointListMortar` (`dynamic_search`) |
| `CheckActivity` | `bool (ModelPart&, const bool ThrowError = true)` | yes | `true` if at least one `SLAVE` node is `ACTIVE`; raises `CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?` when `ThrowError` and none is | contact solvers (`ensure_contact`), `ExplicitPenaltyContactProcess` (`CONTACT` flag) |
| `CheckModelPartHasRotationDoF` | `bool (ModelPart&)` | no | `true` if the elements carry `ROTATION` DoFs | displacement contact criteria (`ROTATION_DOF_IS_CONSIDERED`) |
| `CleanContactModelParts` | `void (ModelPart&)` | yes | Flags `TO_ERASE` every condition whose geometry has geometry parts (i.e. paired conditions) and removes them from all levels | `SearchBaseProcess` (remeshing), `ContactRemeshMmgProcess` |
| `ComputeExplicitContributionConditions` | `void (ModelPart&)` | yes | Calls `AddExplicitContribution(ProcessInfo)` on every condition of the model part (typically `ComputingContact`), which integrates `WEIGHTED_GAP`, `WEIGHTED_SLIP`, `NODAL_AREA` | search, criteria (`PreCriteria`/`PostCriteria`), strategies, `ExplicitPenaltyContactProcess` |
| `ActivateConditionWithActiveNodes` | `void (ModelPart&)` | yes | Sets a paired condition `ACTIVE` if any node of its slave (parent) geometry is `ACTIVE` | `ExplicitPenaltyContactProcess` |
| `GetHalfJumpCenter` | `array_1d<double,3> (GeometryType&)` | no | Center of the geometry displaced by half the Newmark jump (uses `DELTA_COORDINATES`) | search point list in dynamic problems |
| `ComputeTangentMatrixSlip<TDim,TNumNodes>` | `BoundedMatrix<double,TNumNodes,TDim> (const GeometryType&, const std::size_t StepSlip = 1)` | no | Nodal tangent matrix built from the direction of `WEIGHTED_SLIP` at the given buffer step | `PenaltyFrictionalMortarContactCondition` |
| `GetVariableMatrix` (private) | `Matrix (const GeometryType&, const Variable<array_1d<double,3>>&)` | no | Nodal values of a vector variable as a matrix | `ComputeStepJump` |

Test: `tests/cpp_tests/utilities/test_contact_utilities.cpp` (`CheckModelPartHasRotationDoF`).

## `ActiveSetUtilities`

[`active_set_utilities.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_utilities/active_set_utilities.h) implements the update of the active set of the semi-smooth Newton method (thesis §4.3.3.5): after each linear solve the convergence criteria evaluate the *augmented* contact pressure of every slave node, activate the nodes in compression, deactivate the nodes in tension and, in frictional problems, switch between stick and slip. The five functions share the same structure:

1. they run only when the model part is flagged `INTERACTION` (full semi-smooth Newton) **or** `NL_ITERATION_NUMBER == 1` (simplified semi-smooth Newton, in which the active set is frozen after the first iteration; see [`INTERACTION` and `simplified_semi_smooth_newton`](Architecture.html#interaction-and-simplified_semi_smooth_newton));
2. they loop in parallel over the nodes of the `Contact` sub-model-part and consider only `SLAVE` nodes;
3. the penalty of the node is the nodal `INITIAL_PENALTY` if present (set by `ALMFastInit`, adapted by `AALMAdaptPenaltyValueProcess` or `ComputeDynamicFactorProcess`), otherwise the global `ProcessInfo[INITIAL_PENALTY]`;
4. the augmented normal pressure $$\bar{\lambda}_n$$ is stored in `AUGMENTED_NORMAL_CONTACT_PRESSURE` (and the tangent one in `AUGMENTED_TANGENT_CONTACT_PRESSURE`) for post-processing;
5. the return value counts the flags that changed; the criteria declare the active set converged when it is zero (`ACTIVE_SET_CONVERGED`).

| Function | Signature | Python | Augmented normal pressure $$\bar{\lambda}_n$$ | Activation test | Return |
|---|---|---|---|---|---|
| `ComputePenaltyFrictionlessActiveSet` | `std::size_t (ModelPart&)` | yes | $$\varepsilon_i \tilde{g}_{n,i}$$ | $$\bar{\lambda}_n \lt 0$$ → `ACTIVE` | number of `ACTIVE` changes |
| `ComputePenaltyFrictionalActiveSet` | `array_1d<std::size_t,2> (ModelPart&, const bool PureSlip = false, const SizeType EchoLevel = 0)` | yes | $$\varepsilon_i \tilde{g}_{n,i}$$ | as above; stick if $$\Vert \bar{\boldsymbol{\lambda}}_t \Vert \le -\mu \bar{\lambda}_n$$ with $$\bar{\boldsymbol{\lambda}}_t = \tau \varepsilon_i \tilde{\mathbf{g}}_{t,i}$$ | `[0]` `ACTIVE` changes, `[1]` `SLIP` changes |
| `ComputeALMFrictionlessActiveSet` | `std::size_t (ModelPart&)` | no | $$k \lambda_{n,i} + \varepsilon_i \tilde{g}_{n,i}$$ with $$\lambda_{n,i}$$ = `LAGRANGE_MULTIPLIER_CONTACT_PRESSURE` | $$\bar{\lambda}_n \lt 0$$ → `ACTIVE`, and on activation $$\lambda_{n,i} \leftarrow \bar{\lambda}_n / k$$ | number of `ACTIVE` changes |
| `ComputeALMFrictionlessComponentsActiveSet` | `std::size_t (ModelPart&)` | no | $$k\, \boldsymbol{\lambda}_i \cdot \mathbf{n}_i + \varepsilon_i \tilde{g}_{n,i}$$ with $$\boldsymbol{\lambda}_i$$ = `VECTOR_LAGRANGE_MULTIPLIER` | as above; on activation $$\boldsymbol{\lambda}_i \leftarrow \mathbf{n}_i \bar{\lambda}_n / k$$ | number of `ACTIVE` changes |
| `ComputeALMFrictionalActiveSet` | `array_1d<std::size_t,2> (ModelPart&, const bool PureSlip = false, const SizeType EchoLevel = 0)` | no | $$k\, \boldsymbol{\lambda}_i \cdot \mathbf{n}_i + \varepsilon_i \tilde{g}_{n,i}$$ | as above; slip if $$\Vert \bar{\boldsymbol{\lambda}}_t \Vert / (-\mu \bar{\lambda}_n) \gt \theta$$ | `[0]` `ACTIVE` changes, `[1]` `SLIP` changes |

Here $$\varepsilon_i$$ is the nodal penalty, $$k$$ = `SCALE_FACTOR`, $$\tau$$ = `TANGENT_FACTOR`, $$\mu$$ the nodal `FRICTION_COEFFICIENT`, $$\tilde{g}_{n,i}$$ = `WEIGHTED_GAP` and $$\tilde{\mathbf{g}}_{t,i}$$ = `WEIGHTED_SLIP`. The comparison $$\bar{\lambda}_n \lt 0$$ is strict (a source comment flags `<` versus `<=` as a possible point of discussion); a node with $$\bar{\lambda}_n \ge 0$$ that was `ACTIVE` is deactivated (frictional variants also zero `WEIGHTED_SLIP` and reset `SLIP`).

**Frictional details.** In `ComputeALMFrictionalActiveSet` the tangent multiplier is $$\boldsymbol{\lambda}_{t} = \boldsymbol{\lambda} - (\boldsymbol{\lambda}\cdot\mathbf{n})\mathbf{n}$$ and the augmented tangent pressure is

<p align="center">$$ \bar{\boldsymbol{\lambda}}_t = k\, \boldsymbol{\lambda}_t + c\, \tau\, \varepsilon_i\, \tilde{\mathbf{g}}_{t,i}, \qquad c = \begin{cases} \text{SLIP\_AUGMENTATION\_COEFFICIENT} & \text{node currently } \texttt{SLIP} \\ 1 & \text{node currently stick} \end{cases} $$</p>

The stick/slip threshold is $$\theta = 1 - \text{SLIP\_THRESHOLD}$$ for a node currently in slip and $$\theta = 1$$ otherwise, which introduces a hysteresis that prevents chattering between the two states (`slip_threshold`, default `2.0e-2`, and `slip_augmentation_coefficient`, default `0.0`, of `ALMContactProcess`). A node found in slip gets `AUGMENTED_TANGENT_CONTACT_PRESSURE` $$= -\mu \bar{\lambda}_n \, \boldsymbol{\lambda}_t / \Vert \boldsymbol{\lambda}_t \Vert$$ (the Coulomb limit along the tangent multiplier direction); in the penalty variant the direction is that of `WEIGHTED_SLIP`. A newly activated node receives $$\boldsymbol{\lambda}_i \leftarrow \mathbf{n}_i \bar{\lambda}_n / k + \bar{\boldsymbol{\lambda}}_t / k$$ (the tangent part only when $$\mu \gt 0$$). With `PureSlip = true` the `SLIP` flag is always set and never counted as a change; nodes that should stick are reported with a warning when `EchoLevel > 0`. The frictional convergence criteria pass their `PURE_SLIP` option and echo level (`ALMFrictionalMortarConvergenceCriteria`, `PenaltyFrictionalMortarConvergenceCriteria`).

**Callers.** `ALMFrictionlessMortarConvergenceCriteria`, `ALMFrictionlessComponentsMortarConvergenceCriteria`, `ALMFrictionalMortarConvergenceCriteria`, `PenaltyFrictionlessMortarConvergenceCriteria` and `PenaltyFrictionalMortarConvergenceCriteria` call the matching function in `PostCriteria`; `ExplicitPenaltyContactProcess.ExecuteInitializeSolutionStep` calls the two penalty functions from Python, which is why only those two are exposed.

Tests: `tests/cpp_tests/utilities/test_active_set_utilities.cpp` — `ComputePenaltyFrictionlessActiveSet`, `ComputePenaltyFrictionalActiveSet`, `ComputeALMFrictionlessActiveSet`, `ComputeALMFrictionlessComponentsActiveSet`, `ComputeALMFrictionalActiveSet` (three-node model parts with prescribed gaps and multipliers; the expected return value and the resulting `ACTIVE`/`SLIP` flags are checked).

## `SelfContactUtilities`

[`self_contact_utilities.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_utilities/self_contact_utilities.h) contains the automatic master/slave assignment used when `assume_master_slave["N"]` is empty (`predefined_master_slave = false`), described in [Self-contact](../Contact_Search/Self_Contact.html).

| Function | Signature | Python | Purpose |
|---|---|---|---|
| `ComputeSelfContactPairing` | `void (ModelPart&, const std::size_t EchoLevel = 0)` | yes | Resets the `MASTER`/`SLAVE`/`ACTIVE` flags of the conditions and nodes, orders the conditions by proximity and walks them: a condition not yet `MASTER` becomes `SLAVE` and every candidate in its `INDEX_MAP` that is not yet defined, does not share nodes with it and whose nodes are not already `SLAVE` becomes `MASTER` (its own `INDEX_MAP` is cleared); node flags follow the conditions. Conditions ending with both `MASTER` and `SLAVE` nodes are forced to `MASTER` with a warning (`EchoLevel > 0`); an inconsistent count raises |
| `FullAssignmentOfPairs` | `void (ModelPart&)` | yes | Brute force: gives every condition an `INDEX_MAP` containing all the other conditions (used to test the pairing without a tree search) |
| `NotPredefinedMasterSlave` | `void (ModelPart&)` | yes | For the current pairing, sets `SLAVE` on the conditions that have entries in their `INDEX_MAP` and `MASTER` on the others, and propagates the flags to the nodes |

`BaseContactSearchProcess::UpdateMortarConditions` calls `NotPredefinedMasterSlave` before the tree search and `ComputeSelfContactPairing` after it when the roles are not predefined. Tests: `tests/cpp_tests/utilities/test_selfcontact_utilities.cpp` (`SelfContactUtilities1`, `SelfContactUtilities2`, `SelfContactUtilities3`: `FullAssignmentOfPairs` followed by `ComputeSelfContactPairing` on 2D and 3D model parts).

## `DerivativesUtilities`

[`derivatives_utilities.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_utilities/derivatives_utilities.h) computes the directional derivatives of the mortar integration with respect to the nodal displacements that the consistently linearized contact conditions need: derivatives of the Jacobian, of the shape functions (through the integration-cell vertices), of the dual shape functions ($$\mathbf{A}_e$$) and of the normals. The theory and the equation-by-equation correspondence are in [Linearisation and derivatives](../Theory/Linearisation_And_Derivatives.html); this section is the API reference.

The class is a static template `DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster = TNumNodes>` explicitly instantiated for the five geometry pairs (`<2,2>`, `<3,3>`, `<3,4>`, `<3,3,4>`, `<3,4,3>`) times `TFrictional` $$\in$$ {`false`, `true`} times `TNormalVariation` $$\in$$ {`false`, `true`}, i.e. 20 instantiations. Its type aliases connect it to the core mortar classes: `DerivativeDataType` is `DerivativeData` or `DerivativeDataFrictional`, `GeneralVariables` is `MortarKinematicVariablesWithDerivatives`, `AeData` is `DualLagrangeMultiplierOperatorsWithDerivatives`, `MortarConditionMatrices` is `MortarOperatorWithDerivatives`, `BelongType` is the `PointBelongs*` enum of the pair and `DecompositionType` is `Line2D2<Point>` or `Triangle3D3<Point>` (the integration cell). Constants: `CheckThresholdCoefficient = 1.0e-12`, `ZeroTolerance = machine epsilon`.

| Method | Purpose |
|---|---|
| `CalculateDeltaDetjSlave(DecompGeom, rVariables, rDerivativeData)` | Derivative of the Jacobian determinant of the integration cell (thesis §4.6.1.1 / §4.6.2.1) |
| `GPDeltaNormalSlave(rJacobian, rDNDe)`, `GPDeltaNormalMaster(rJacobian, rDNDe)` | Local (Gauss-point) increment of the normal on each side |
| `DeltaNormalCenter(rThisGeometry)` | Increment of the normal at the center of a geometry |
| `CalculateDeltaNormalSlave(rDeltaNormal, rThisGeometry)`, `CalculateDeltaNormalMaster(rDeltaNormal, rThisGeometry)` | Nodal normal derivatives on slave and master (`TNormalVariation`) (§4.6.1.4 / §4.6.2.4) |
| `CalculateDeltaCellVertex(rVariables, rDerivativeData, rTheseBelongs, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, rNormal)` | Derivatives of the vertices of the integration cell; uses `LocalDeltaVertex`, `LocalDeltaSegmentN1`, `DeltaPointLocalCoordinatesSlave/Master` and `ConvertAuxHashIndex` (the `PointBelongs` hash tells to which cut each vertex belongs) |
| `CalculateDeltaN1(...)`, `CalculateDeltaN(...)` | Derivatives of the slave (`N1`) and of both slave and master (`N`) shape functions (§4.6.1.2 / §4.6.2.2) |
| `CalculateDeltaPosition(...)` (four overloads) | Increment of displacements between the current and previous iteration as a matrix, a matrix for a given DoF, a vector for a node or a scalar for a node/component |
| `CalculateAeAndDeltaAe(rSlaveGeometry, rSlaveNormal, rMasterGeometry, rDerivativeData, rVariables, ConsiderNormalVariation, rIntegrationUtility, rConditionsPointsSlave, ThisIntegrationMethod, AxiSymCoeff)` | Computes the dual operator $$\mathbf{A}_e$$ and its derivatives $$\Delta\mathbf{A}_e$$ over the segmented cells; returns `false` when the integration area is null (§4.6.1.3 / §4.6.2.3) |
| `LocalDeltaVertex(...)` | Derivative of one cell vertex in local terms |
| `ComputeRenormalizerMatrix(rDiffVector, rDeltaNormal)` | Auxiliary matrix that keeps the normal unitary after an increment |
| `PreviousNormalGeometry(rThisGeometry, rDeltaNormal)` | Normal in the previous configuration |

The helper class `ImplementationDerivativesUtilities` (same header) holds the non-templated implementation of `ComputeRenormalizerMatrix`. The `MortarContactCondition` family calls these methods from `CalculateLocalLHS`/`CalculateConditionSystem` when `TNormalVariation` is `true` or the LHS is assembled; the automatically generated conditions receive the derivatives through `DerivativeData` (`DeltaDetjSlave`, `DeltaPhi`, `DeltaN1`, `DeltaN2`, `DeltaNormalSlave`, `DeltaNormalMaster`, `DeltaCellVertex`).

Tests: `tests/cpp_tests/utilities/test_derivatives_utilities.cpp` contains 49 cases that compare the analytical derivatives with finite differences for `DerivativesUtilities<TDim, TNumNodes, false, true>`: `JacobianDerivatives{Line1-3, Triangle1-6, Quadrilateral1-3}`, `ShapeFunctionDerivatives{Line1-4, Triangle1-6, Quadrilateral1-3}`, `DualShapeFunctionDerivatives{Line1-3, Triangle1-6, Quadrilateral1-3}` and `NormalDerivatives{Line1-3, Triangle1-6, Quadrilateral1-3}` (`CalculateDeltaNormalSlave/Master`, `CalculateAeAndDeltaAe`, `CalculateDeltaCellVertex`, `CalculateDeltaDetjSlave`, `CalculateDeltaN`, `DeltaNormalCenter`).

## `MortarExplicitContributionUtilities` and `AuxiliaryOperationsUtilities`

[`mortar_explicit_contribution_utilities.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_utilities/mortar_explicit_contribution_utilities.h) evaluates the mortar operators of a paired condition *without* derivatives, i.e. what is needed for the residual-only ("explicit") contribution that feeds the weighted gap, the weighted slip and the nodal area used by the search and by the active-set update. The static template `MortarExplicitContributionUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster = TNumNodes>` takes a `FrictionalCase` enum value (`FRICTIONLESS = 0`, `FRICTIONLESS_COMPONENTS = 1`, `FRICTIONAL = 2`, `FRICTIONLESS_PENALTY = 3`, `FRICTIONAL_PENALTY = 4`) as `TFrictional`; it is instantiated 50 times (5 pairs × 5 cases × 2 normal-variation values). Two compile-time constants mirror the conditions: `MatrixSize` (the size of the local system: $$d(n_m + n_s) + n_s$$ for `FRICTIONLESS`, $$d(n_m + 2 n_s)$$ for `FRICTIONLESS_COMPONENTS` and `FRICTIONAL`, $$d(n_m + n_s)$$ for the penalty cases) and `IsFrictional` (`true` for `FRICTIONAL` and `FRICTIONAL_PENALTY`).

| Method | Purpose |
|---|---|
| `AddExplicitContributionOfMortarCondition(pCondition, rCurrentProcessInfo, IntegrationOrder = 2, AxisymmetricCase = false, ComputeNodalArea = false, ComputeDualLM = true, rAreaVariable = NODAL_AREA)` | Integrates the mortar operators $$\mathbf{D}$$, $$\mathbf{M}$$ of the pair with `ExactMortarIntegrationUtility`, adds $$\sum_j (D_{ij} x^s_j - M_{ij} x^m_j)\cdot \mathbf{n}_i$$ to `WEIGHTED_GAP` of the slave nodes (atomic adds) and, optionally, $$D_{ii}$$ to the area variable; returns the operators |
| `AddExplicitContributionOfMortarFrictionalCondition(pCondition, rCurrentProcessInfo, rPreviousMortarOperators, IntegrationOrder = 2, AxisymmetricCase = false, ComputeNodalArea = false, ComputeDualLM = true, rAreaVariable = NODAL_AREA, ConsiderObjetiveFormulation = false)` | As above plus the weighted slip `WEIGHTED_SLIP`, computed from the current and previous mortar operators (objective formulation optional) |
| `ExplicitCalculateAe(rSlaveGeometry, rVariables, rConditionsPointsSlave, rAe, rIntegrationMethod, AxiSymCoeff = 1.0)` | Dual-multiplier operator $$\mathbf{A}_e$$ without derivatives |
| `ExplicitCalculateKinematics(pCondition, rVariables, rAe, rNormalMaster, rLocalPointDecomp, rLocalPointParent, rGeometryDecomp, DualLM = true)` | Shape functions, Jacobian and dual shape functions at a cell point, without derivatives |
| `ComputeNodalArea(pCondition, rCurrentProcessInfo, rAreaVariable = NODAL_AREA, IntegrationOrder = 2, AxisymmetricCase = false)` | Adds the diagonal of $$\mathbf{D}$$ to the nodal area variable |
| `ComputePreviousMortarOperators(pCondition, rCurrentProcessInfo, rPreviousMortarOperators, IntegrationOrder = 2, AxisymmetricCase = false, ComputeDualLM = true)` | Mortar operators in the previous configuration (frictional slip) |
| `CalculateKinematics(pCondition, rVariables, rDerivativeData, rNormalMaster, rLocalPointDecomp, rLocalPointParent, rGeometryDecomp, DualLM = true)` | Kinematics with the `DerivativeData` container (shared with the implicit path) |
| `MasterShapeFunctionValue(pCondition, rVariables, rNormalMaster, rLocalPoint)` | Projects a slave point on the master and evaluates the master shape functions |

The namespace `AuxiliaryOperationsUtilities` (same header) provides `GetAxisymmetricCoefficient(pCondition, rNSlave)` (the $$2\pi r$$ weight of the axisymmetric conditions) and `CalculateRadius(pCondition, rNSlave)`. `MortarContactCondition::AddExplicitContribution` and the frictional/penalty conditions delegate to these utilities, and `ContactUtilities::ComputeExplicitContributionConditions` is the model-part-level loop that triggers them. The weighted-gap tests `WeightedGap1`–`WeightedGap9` of `tests/cpp_tests/processes/test_weighted_gap.cpp` exercise this path indirectly (see [Gap computation](../Contact_Search/Gap_Computation.html#weighted-gap-unit-tests-as-executable-definitions)).

## `InterfacePreprocessCondition`

[`interface_preprocess.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_utilities/interface_preprocess.h) turns a user model part that marks the contact interface into a model part of surface/line **conditions** ready for the search. It is constructed with the main model part (`InterfacePreprocessCondition(ModelPart& rMainModelPart)`) and has one public method, `GenerateInterfacePart(ModelPart& rInterfacePart, Parameters ThisParameters)`, both exposed to Python. Defaults:

```json
{
    "simplify_geometry"                    : false,
    "contact_property_id"                  : 0
}
```

`GenerateInterfacePart` distinguishes three situations:

1. **The interface model part already has conditions** (recommended input). Nothing is created; `CheckAndCreateProperties` verifies that the properties of the conditions carry `YOUNG_MODULUS`. If they do not, a new `Properties` (id = number of properties + 1) is created, the parent element of the first condition is located by comparing sorted node ids with the element faces (`CheckOnTheFace`, `GenerateBoundariesEntities`), and `YOUNG_MODULUS`, `THICKNESS` (optional) and, in frictional problems (main model part `SLIP`), `FRICTION_COEFFICIENT` are copied from the element properties (`CopyProperties`) into the new properties, which are then assigned to every interface condition.
2. **The interface model part has only nodes.** The conditions of the main model part are renumbered consecutively (`ReorderConditions`, returns the next free id). If `contact_property_id == 0`, `CreateNewProperties` creates one new `Properties` per existing `Properties` (ids continue after the existing ones) copying `YOUNG_MODULUS`, `THICKNESS` and, with a deprecation warning, `FRICTION_COEFFICIENT`; otherwise the single property `contact_property_id` is created in the main model part. Then every element is visited: for solid elements (`WorkingSpaceDimension == LocalSpaceDimension`) each boundary entity is passed to `GenerateEdgeCondition` (2D) or `GenerateFaceCondition` (3D); for structural elements (shells, beams) the element geometry itself is passed. A boundary generates a condition only when **all** its nodes are `INTERFACE`. The condition type is `LineCondition2D2N` (or `LineCondition2D3N` for three-node edges) in 2D and `SurfaceCondition3D3N` / `3D4N` / `3D6N` / `3D8N` / `3D9N` in 3D according to the number of nodes; `simplify_geometry: true` forces `LineCondition2D2N` / `SurfaceCondition3D3N` (linear geometries). `CreateNewCondition` creates the condition in the interface model part and `AssignMasterSlaveCondition` sets it `SLAVE` if all its nodes are `SLAVE` and `MASTER` otherwise.
3. **Neither conditions nor nodes**: `KRATOS_ERROR` ("Nor conditions or nodes on the interface. Check your flags").

Finally `PrintNodesAndConditions` reports the number of interface nodes and conditions and raises if either is zero. `SearchBaseProcess._interface_preprocess` calls it for every user model part of the pair with `contact_property_id = search_property_ids[N]` (the ALM/penalty/MPC processes name it `contact_property_ids`, the mesh-tying process `mesh_tying_property_ids`), after which the conditions are transferred to `ContactSub<N>`. Tests: `tests/cpp_tests/utilities/test_interface_preprocess_utilities.cpp` (`InterfacePreprocessCondition2D`, `InterfacePreprocessCondition3D`); `tests/test_dynamic_search.py` also uses it from Python.

## `ProcessFactoryUtility`

[`process_factory_utility.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_python/process_factory_utility.h) lives in `custom_python/` because it depends on pybind11: it holds a `std::vector<pybind11::object>` of *Python* processes and lets the C++ strategies call their hooks. Constructors: default, `ProcessFactoryUtility(py::list& ProcessesList)` and `ProcessFactoryUtility(py::object& rProcess)`.

| Method | Python | Behavior |
|---|---|---|
| `AddProcess(rProcess)`, `AddProcesses(ProcessesList)` | yes | Append one process or a list |
| `ExecuteMethod(rNameMethod)` | yes | For every stored object, `process.attr(rNameMethod)()` — any zero-argument method name can be invoked |
| `ExecuteInitialize`, `ExecuteBeforeSolutionLoop`, `ExecuteInitializeSolutionStep`, `ExecuteFinalizeSolutionStep`, `ExecuteBeforeOutputStep`, `ExecuteAfterOutputStep`, `ExecuteFinalize`, `IsOutputStep`, `PrintOutput`, `Clear` | yes | `ExecuteMethod("<name>")` |
| `Info`, `PrintInfo`, `PrintData` | no | Print the number of stored processes |

The contact solvers (`contact_structural_mechanics_static_solver.py`, `contact_structural_mechanics_implicit_dynamic_solver.py`) implement `AddProcessesList(processes_list)` and `AddPostProcess(post_process)`, wrapping the Python lists in `CSMA.ProcessFactoryUtility`, and pass the two objects to `ResidualBasedNewtonRaphsonContactStrategy` / `LineSearchContactStrategy` through `auxiliary_methods_solvers.AuxiliaryNewton` / `AuxiliaryLineSearch`. The strategies keep them as `mpMyProcesses` and `mpPostProcesses` and, when `adaptative_strategy` is enabled and a step is split (`AdaptativeStep`), re-run `ExecuteInitializeSolutionStep`, `ExecuteFinalizeSolutionStep`, `ExecuteBeforeOutputStep`, `PrintOutput` and `ExecuteAfterOutputStep` on the sub-steps, which is not possible from the analysis stage alone. A warning is printed when the strategy has no process list and `echo_level > 0`. Tests: `tests/test_process_factory.py` (`test_process_factory`, `test_processes_list_factory`).

## `logging_settings.hpp`

[`logging_settings.hpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ContactStructuralMechanicsApplication/custom_utilities/logging_settings.hpp) is a header of preprocessor macros for colored console debugging, included by `mortar_contact_condition.h`, `mesh_tying_mortar_condition.h` and `mixedulm_linear_solver.h`. It defines ANSI color codes (`RESET`, `RED`, `GREEN`, `YELLOW`, `BLUE`, `MAGENTA`, `CYAN`, `WHITE`, `DK_GREY`, `LT_*`, `BOLD`, `DIM`, `UNDERLINE`), generic message macros (`LOG_GENERAL`, `LOG_DEBUG`, `LOG_INFO`, `LOG_WARNING`, `*_HEADER` variants, `DEBUG_MSG`, `INFO_MSG`, `WARNING_MSG`, `ERROR_MSG`), pretty printers for tensors (`LOG_MATRIX_PRETTY`, `LOG_VECTOR_PRETTY`, `LOG_VECTOR3`, `LOG_VECTOR2`, `LOG_SCALAR`) and a condition dump header (`LOG_CONDITION_HEADER(master, slave)`), configured by `RESET_LOG_SETTINGS`, `TENSOR_LOG_SETTINGS` and `CONDITION_LOG_SETTINGS`. They are debugging aids with no effect on the simulation and no Python exposure.

## Kratos-core utilities the application depends on

The mortar machinery itself (segmentation, integration, operators, dual shape functions) lives in the Kratos core so that the mapper (`SimpleMortarMapperProcess`) and this application share it. The application uses three headers.

### `MortarUtilities` — [`kratos/utilities/mortar_utilities.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/utilities/mortar_utilities.h)

| Function | Used in |
|---|---|
| `LengthCheck(rGeometryLine, Tolerance)`, `HeronCheck(rGeometryTriangle)` | Conditions and search: discard degenerate integration cells (zero length / Heron area) |
| `RotatePoint(rPointToRotate, rPointReferenceRotation, rSlaveTangentXi, rSlaveTangentEta, Inversed)` | Projection of the master onto the slave tangent plane (tested in `TestCheckRotation`) |
| `GaussPointUnitNormal(rN, rGeometry)` | Normal at a Gauss point from the nodal normals (conditions) |
| `ComputeNodesMeanNormalModelPart(rModelPart, ComputeConditions)` | Nodal normals of the contact model part (search, criteria) |
| `ComputeNodesTangentModelPart(rModelPart, pSlipVariable, SlipCoefficient, SlipAlways)` | Nodal tangents from the multiplier or from `WEIGHTED_SLIP` (`BaseMortarConvergenceCriteria::PreCriteria`, `ALMContactProcess._reset_slip_flag`) |
| `InvertNormalForFlag<TContainerType>(rContainer, rFlag)` | `NormalCheckProcess` (inverts flagged elements/conditions) |
| `GetCoordinates<TDim,TNumNodes>(rGeometry, Current, Step)` | Nodal coordinates as a matrix (conditions, derivatives) |
| `ComputeTangentMatrix<TDim,TNumNodes>(rGeometry)` | Nodal tangent matrix (frictional conditions) |
| `GetVariableVector<TNumNodes>(rGeometry, rVariable, Step)`, `GetVariableMatrix<TDim,TNumNodes>(rGeometry, rVariable, Step)` | Nodal values of scalar/vector variables (`WEIGHTED_GAP`, `LAGRANGE_MULTIPLIER_CONTACT_PRESSURE`, `VECTOR_LAGRANGE_MULTIPLIER`, `NORMAL`, …); the most frequently used helpers in the generated condition code |

### `ExactMortarIntegrationUtility` — [`kratos/utilities/exact_mortar_segmentation_utility.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/utilities/exact_mortar_segmentation_utility.h)

`ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>` performs the exact segmentation of a slave geometry against a master geometry (clipping in 2D, polygon clipping and triangulation in 3D) and returns the integration cells as `ConditionArrayListType` (arrays of `PointBelong` points when `TBelong` is `true`, which the derivatives need). Constructor: `(IntegrationOrder = 0, DistanceThreshold = max, EchoLevel = 0, ZeroToleranceFactor = 1.0, ConsiderDelaunator = false)`. The application uses `GetExactIntegration(rOriginalSlaveGeometry, rSlaveNormal, rOriginalMasterGeometry, rMasterNormal, rConditionsPointsSlave)` in the conditions and the explicit utilities, `GetTotalArea` and `TestGetExactAreaIntegration` for debugging (`SearchBaseProcess.__get_integration_area` calls the Python classes `KM.ExactMortarIntegrationUtility2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N`), and `SetEchoLevel`. `INTEGRATION_ORDER_CONTACT`, `DISTANCE_THRESHOLD`, `ZERO_TOLERANCE_FACTOR` and `CONSIDER_TESSELLATION` are the application-side knobs forwarded to it. The theory is in [Mortar integration and dual Lagrange multipliers](../Theory/Mortar_Integration_And_Dual_Lagrange_Multipliers.html); `tests/cpp_tests/utilities/test_integration_utilities.cpp` (`MassMatrixIntegrationTriangle`, `MassMatrixIntegrationQuadrilateral`, `MassMatrixIntegrationQuadrilateralDeformed`, `TestCheckRotation`) verifies it from the application side.

### `mortar_classes.h` — [`kratos/includes/mortar_classes.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/includes/mortar_classes.h)

| Type | Content | Application use |
|---|---|---|
| `PointBelongs`, `PointBelongsLine2D2N`, `PointBelongsTriangle3D3N`, `PointBelongsQuadrilateral3D4N`, `PointBelongsTriangle3D3NQuadrilateral3D4N`, `PointBelongsQuadrilateral3D4NTriangle3D3N` | Enum "hashes" telling to which slave node, master node or edge intersection a cell vertex belongs | `DerivativesUtilities::BelongType`, `CalculateDeltaCellVertex` |
| `PointBelong<TNumNodes, TNumNodesMaster>` | A `Point` with a `PointBelongs*` tag | Integration cells with derivatives |
| `MortarKinematicVariables<TNumNodes, TNumNodesMaster>` | $$N$$, $$\Phi$$, Jacobian, integration weight at a cell point | Conditions, explicit utilities |
| `MortarKinematicVariablesWithDerivatives<TDim, TNumNodes, TNumNodesMaster>` | Adds the shape-function gradients and Jacobian matrices | `DerivativesUtilities::GeneralVariables` |
| `DerivativeData<TDim, TNumNodes, TNumNodesMaster>`, `DerivativeDataFrictional<…>` | Containers for $$\Delta\det J$$, $$\Delta N_1$$, $$\Delta N_2$$, $$\Delta\Phi$$, $$\Delta\mathbf{n}$$, $$\Delta$$ cell vertices, penalty, scale factor, normals; the frictional one adds the previous-step data | `DerivativesUtilities::DerivativeDataType`, condition `CalculateLocalLHS/RHS` |
| `MortarOperator<TNumNodes, TNumNodesMaster>` | $$\mathbf{D}$$ and $$\mathbf{M}$$ (Popp's definition) with `CalculateMortarOperators` and `ComputePOperator` | Explicit utilities, mesh tying, MPC condition |
| `MortarOperatorWithDerivatives<TDim, TNumNodes, TFrictional, TNumNodesMaster>` | Adds $$\Delta\mathbf{D}$$, $$\Delta\mathbf{M}$$ | `DerivativesUtilities::MortarConditionMatrices` |
| `DualLagrangeMultiplierOperators<TNumNodes, TNumNodesMaster>` | $$\mathbf{M}_e$$, $$\mathbf{D}_e$$ and $$\mathbf{A}_e = \mathbf{D}_e \mathbf{M}_e^{-1}$$ of the dual basis | `ExplicitCalculateAe` |
| `DualLagrangeMultiplierOperatorsWithDerivatives<TDim, TNumNodes, TFrictional, TNumNodesMaster>` | Adds $$\Delta\mathbf{A}_e$$ | `DerivativesUtilities::AeData`, `CalculateAeAndDeltaAe` |

## Utilities → tests

All C++ tests belong to the suite `KratosContactStructuralMechanicsFastSuite` (`tests/cpp_tests/contact_structural_mechanics_fast_suite.h`) and are run with the `KratosContactStructuralMechanicsCoreTest` executable or `python run_cpp_unit_tests.py`.

| Utility | Test file (`tests/cpp_tests/utilities/`) | Cases |
|---|---|---|
| `ActiveSetUtilities` | `test_active_set_utilities.cpp` | `ComputePenaltyFrictionlessActiveSet`, `ComputePenaltyFrictionalActiveSet`, `ComputeALMFrictionlessActiveSet`, `ComputeALMFrictionlessComponentsActiveSet`, `ComputeALMFrictionalActiveSet` |
| `ContactUtilities` | `test_contact_utilities.cpp` | `CheckModelPartHasRotationDoF` |
| `DerivativesUtilities` | `test_derivatives_utilities.cpp` | 49 cases: `DualShapeFunctionDerivativesLine1-3`, `JacobianDerivativesLine1-3`, `NormalDerivativesLine1-3`, `ShapeFunctionDerivativesLine1-4`, `JacobianDerivativesTriangle1-6`, `ShapeFunctionDerivativesTriangle1-6`, `DualShapeFunctionDerivativesTriangle1-6`, `NormalDerivativesTriangle1-6`, `JacobianDerivativesQuadrilateral1-3`, `ShapeFunctionDerivativesQuadrilateral1-3`, `DualShapeFunctionDerivativesQuadrilateral1-3`, `NormalDerivativesQuadrilateral1-3` |
| `ExactMortarIntegrationUtility`, `MortarUtilities::RotatePoint` (core) | `test_integration_utilities.cpp` | `MassMatrixIntegrationTriangle`, `MassMatrixIntegrationQuadrilateral`, `MassMatrixIntegrationQuadrilateralDeformed`, `TestCheckRotation` |
| `InterfacePreprocessCondition` | `test_interface_preprocess_utilities.cpp` | `InterfacePreprocessCondition2D`, `InterfacePreprocessCondition3D` |
| `SelfContactUtilities` | `test_selfcontact_utilities.cpp` | `SelfContactUtilities1`, `SelfContactUtilities2`, `SelfContactUtilities3` |
| `MortarExplicitContributionUtilities` (indirect) | `../processes/test_weighted_gap.cpp` | `WeightedGap1`–`WeightedGap9` (with `3b`, `4b`) |
| `ProcessFactoryUtility` | `tests/test_process_factory.py` (Python) | `test_process_factory`, `test_processes_list_factory` |

## Notes and limitations

- `ComputeALMFrictionlessActiveSet`, `ComputeALMFrictionlessComponentsActiveSet` and `ComputeALMFrictionalActiveSet` are not exposed to Python; only the convergence criteria call them. Driving an ALM active-set update from a Python script requires the criteria (or a custom binding).
- The active-set functions read the nodal `FRICTION_COEFFICIENT` written by `ALMFastInit` (an average over the attached conditions), not the `FrictionalLaw` classes described in [Frictional laws and MPC constraint](Frictional_Laws_And_MPC_Constraint.html).
- `InterfacePreprocessCondition` only generates conditions from element boundaries whose nodes are **all** `INTERFACE`; interfaces described by nodes on curved or coarse meshes may therefore miss faces, which is why supplying conditions is the recommended input.
- `CleanContactModelParts` identifies paired conditions by the presence of geometry parts (`NumberOfGeometryParts() > 0`); any other condition built on a composite geometry in the same model part would be removed as well.
- `DerivativesUtilities` and `MortarExplicitContributionUtilities` are compiled for the five geometry pairs only; higher-order interface geometries (`LineCondition2D3N`, `SurfaceCondition3D6N`, …, which `InterfacePreprocessCondition` can create) are not supported by the mortar conditions and must be avoided with `simplify_geometry: true` or by meshing the interface with linear elements.
