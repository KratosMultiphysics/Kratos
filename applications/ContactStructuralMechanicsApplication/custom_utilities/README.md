# Utilities

Compiled helper classes and namespaces shared by the conditions, processes and criteria of the application.

| File | Class / namespace | Purpose | Exposed to Python |
|---|---|---|---|
| `contact_utilities.h/.cpp` | `ContactUtilities` (static class) | Mesh-size helpers (`CalculateRelativeSizeMesh`, `CalculateMaxNodalH`, `CalculateMeanNodalH`, `CalculateMinimalNodalH`), `ComputeExplicitContributionConditions` (residual-only evaluation of the pairs, used by `Predict()` and by the criteria to refresh the weighted gap), `CheckActivity`, `CleanContactModelParts`, `ActivateConditionWithActiveNodes`, `CheckModelPartHasRotationDoF`, `ComputeTangentMatrixSlip`, `GetHalfJumpCenter`. | yes (subset) |
| `active_set_utilities.h/.cpp` | `ActiveSetUtilities` (namespace) | The semi-smooth Newton set updates: `ComputePenaltyFrictionlessActiveSet`, `ComputePenaltyFrictionalActiveSet`, `ComputeALMFrictionlessActiveSet`, `ComputeALMFrictionlessComponentsActiveSet`, `ComputeALMFrictionalActiveSet`. A node is `ACTIVE` when the augmented pressure $\bar\lambda_n = k\lambda_n + \varepsilon\tilde g_n$ (or $\varepsilon\tilde g_n$ for penalty) is negative; frictional variants also flip the `SLIP` flag with the Coulomb check and return the number of active-set and stick/slip changes. | penalty variants |
| `self_contact_utilities.h/.cpp` | `SelfContactUtilities` (namespace) | Automatic master/slave assignment for self-contact (thesis Algorithm 4): `ComputeSelfContactPairing`, `FullAssignmentOfPairs`, `NotPredefinedMasterSlave`. | yes |
| `derivatives_utilities.h/.cpp` | `DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>` | Directional derivatives needed for the consistent tangent of the mortar conditions (thesis §4.6): Jacobian (`CalculateDeltaDetjSlave`), normals (`CalculateDeltaNormalSlave/Master`, `GPDeltaNormalSlave/Master`, `DeltaNormalCenter`), clipping vertices (`CalculateDeltaCellVertex`), shape functions (`CalculateDeltaN1`, `CalculateDeltaN`), positions (`CalculateDeltaPosition`), dual coefficients (`CalculateAeAndDeltaAe`). These are the "AD exceptions" referenced by the generated condition code. | no |
| `mortar_explicit_contribution_utilities.h/.cpp` | `MortarExplicitContributionUtilities<…>`, `AuxiliaryOperationsUtilities` | Residual-only (explicit) mortar contributions: `AddExplicitContributionOfMortarCondition`, `AddExplicitContributionOfMortarFrictionalCondition`, `ExplicitCalculateAe`, `ExplicitCalculateKinematics`, `ComputeNodalArea`, `ComputePreviousMortarOperators`, `CalculateKinematics`, `MasterShapeFunctionValue`; axisymmetric helpers `GetAxisymmetricCoefficient`, `CalculateRadius`. | no |
| `interface_preprocess.h/.cpp` | `InterfacePreprocessCondition` | Creates the interface (skin) conditions of a body from its elements (`GenerateInterfacePart`; parameters `simplify_geometry`, `contact_property_id`). | yes |
| `logging_settings.hpp` | macros | Coloured debug printing of matrices, vectors and conditions (`KRATOS_WATCH_*`). | – |

The `ProcessFactoryUtility` that lets the C++ strategies drive Python process lists lives in `custom_python/process_factory_utility.h`.

Core utilities the application depends on (in `kratos/`): `MortarUtilities` (`utilities/mortar_utilities.h`: nodal normals and tangents, variable matrices), `ExactMortarIntegrationUtility` (`utilities/exact_mortar_segmentation_utility.h`: segmentation and integration cells), the mortar containers of `includes/mortar_classes.h` (`MortarOperator`, `MortarOperatorWithDerivatives`, `DualLagrangeMultiplierOperators`, `DerivativeData`, `MortarKinematicVariables`) and `SimpleMortarMapperProcess`.

Tests: `tests/cpp_tests/utilities/` — `test_derivatives_utilities.cpp` (49 cases), `test_active_set_utilities.cpp` (5), `test_integration_utilities.cpp` (4), `test_selfcontact_utilities.cpp` (3), `test_interface_preprocess_utilities.cpp` (2), `test_contact_utilities.cpp` (1).

## Full documentation

- [Utilities (implementation reference)](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Utilities.html) · [source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Utilities.md)
- [Linearisation and derivatives](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Linearisation_And_Derivatives.md) · [Self-contact](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/Self_Contact.md) · [Gap computation](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Contact_Search/Gap_Computation.md)
