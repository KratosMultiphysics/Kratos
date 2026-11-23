---
type: C++ Function
title: ThisThread
resource: utilities/openmp_utils.h#L108-L115
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/Allocate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/Deallocate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/AddChunk
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/SetMaximumAndMinimumDistanceValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CheckAndResizeThreadStorage
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityResidualLocalContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityGradientContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityFirstDerivativeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntitySecondDerivativeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculatePreviousTimeStepContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityTimeSchemeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityAuxiliaryVariableContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GatherVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/CalculateRHSContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_steady_scheme/ResidualBasedAdjointSteadyScheme/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_scheme/TemplateAddDynamicsToRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/AddDynamicsToRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_implicit_time_scheme/TCalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_implicit_time_scheme/TCalculateRHSContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/AddDynamicsToRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/elimination_build_dof_array_utility/EliminationBuildDofArrayUtility/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline int ThisThread()`

# Called by

- [Allocate](../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/Allocate.md)
- [Deallocate](../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/Deallocate.md)
- [AddChunk](../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/AddChunk.md)
- [SetMaximumAndMinimumDistanceValues](../../../functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/SetMaximumAndMinimumDistanceValues.md)
- [CalculateSystemContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSystemContributions.md)
- [CheckAndResizeThreadStorage](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CheckAndResizeThreadStorage.md)
- [CalculateEntityResidualLocalContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityResidualLocalContributions.md)
- [CalculateEntityGradientContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityGradientContributions.md)
- [CalculateEntityFirstDerivativeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityFirstDerivativeContributions.md)
- [CalculateEntitySecondDerivativeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntitySecondDerivativeContributions.md)
- [CalculatePreviousTimeStepContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculatePreviousTimeStepContributions.md)
- [CalculateEntityTimeSchemeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions.md)
- [CalculateEntityAuxiliaryVariableContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions.md)
- [UpdateEntityTimeSchemeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityTimeSchemeContributions.md)
- [UpdateEntityAuxiliaryVariableContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityAuxiliaryVariableContributions.md)
- [GatherVariables](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GatherVariables.md)
- [CalculateSystemContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/CalculateSystemContributions.md)
- [CalculateRHSContribution](../../../functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/CalculateRHSContribution.md)
- [CalculateSystemContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_steady_scheme/ResidualBasedAdjointSteadyScheme/CalculateSystemContributions.md)
- [TemplateAddDynamicsToRHS](../../../functions/solving_strategies/schemes/residual_based_bdf_scheme/TemplateAddDynamicsToRHS.md)
- [AddDynamicsToRHS](../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/AddDynamicsToRHS.md)
- [TCalculateSystemContributions](../../../functions/solving_strategies/schemes/residual_based_implicit_time_scheme/TCalculateSystemContributions.md)
- [TCalculateRHSContribution](../../../functions/solving_strategies/schemes/residual_based_implicit_time_scheme/TCalculateRHSContribution.md)
- [AddDynamicsToRHS](../../../functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/AddDynamicsToRHS.md)
- [CalculateLocalSensitivity](../../../functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivity.md)
- [ReadNodesBlock](../../../functions/sources/model_part_io/ModelPartIO/ReadNodesBlock.md)
- [SetUpDofArray](../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpDofArray.md)
- [SetUpDofArray](../../../functions/utilities/dof_utilities/elimination_build_dof_array_utility/EliminationBuildDofArrayUtility/SetUpDofArray.md)