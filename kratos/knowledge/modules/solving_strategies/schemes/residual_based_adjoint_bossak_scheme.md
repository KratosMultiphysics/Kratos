---
type: C++ Module
title: residual_based_adjoint_bossak_scheme
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L1-L1363
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/vector
    resolved_by: tree-sitter
    confidence: exact
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/unordered-set
    resolved_by: tree-sitter
    confidence: exact
  - target: external/functional
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-checks-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-parameters-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-schemes-scheme-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/response-functions-adjoint-response-function-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-indirect-scalar-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-adjoint-extensions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-atomic-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ResidualBasedAdjointBossakScheme](../../../classes/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme.md)
- [ResidualBasedAdjointBossakScheme](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/ResidualBasedAdjointBossakScheme.md)
- [Check](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Check.md)
- [Initialize](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Initialize.md)
- [Predict](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Predict.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/InitializeSolutionStep.md)
- [FinalizeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/FinalizeSolutionStep.md)
- [Update](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Update.md)
- [CalculateSystemContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSystemContributions.md)
- [CalculateLHSContribution](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateLHSContribution.md)
- [CalculateSystemContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSystemContributions-2.md)
- [CalculateLHSContribution](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateLHSContribution-2.md)
- [Clear](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Clear.md)
- [Info](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Info.md)
- [PrintInfo](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/PrintInfo.md)
- [PrintData](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/PrintData.md)
- [BossakConstants](../../../classes/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/BossakConstants.md)
- [CalculateGradientContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateGradientContributions.md)
- [CalculateGradientContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateGradientContributions-2.md)
- [CalculateFirstDerivativeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateFirstDerivativeContributions.md)
- [CalculateFirstDerivativeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateFirstDerivativeContributions-2.md)
- [CalculateSecondDerivativeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSecondDerivativeContributions.md)
- [CalculateSecondDerivativeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSecondDerivativeContributions-2.md)
- [CalculateResidualLocalContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateResidualLocalContributions.md)
- [CalculateResidualLocalContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateResidualLocalContributions-2.md)
- [CalculateTimeSchemeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateTimeSchemeContributions.md)
- [CalculateTimeSchemeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateTimeSchemeContributions-2.md)
- [CalculateAuxiliaryVariableContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateAuxiliaryVariableContributions.md)
- [CalculateAuxiliaryVariableContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateAuxiliaryVariableContributions-2.md)
- [CheckAndResizeThreadStorage](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CheckAndResizeThreadStorage.md)
- [CalculateEntityResidualLocalContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityResidualLocalContributions.md)
- [CalculateEntityGradientContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityGradientContributions.md)
- [CalculateEntityFirstDerivativeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityFirstDerivativeContributions.md)
- [CalculateEntitySecondDerivativeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntitySecondDerivativeContributions.md)
- [CalculatePreviousTimeStepContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculatePreviousTimeStepContributions.md)
- [CalculateEntityTimeSchemeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions.md)
- [CalculateEntityAuxiliaryVariableContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions.md)
- [CalculateNodeNeighbourCount](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateNodeNeighbourCount.md)
- [UpdateTimeSchemeAdjoints](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints.md)
- [UpdateEntityTimeSchemeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityTimeSchemeContributions.md)
- [UpdateTimeSchemeVariablesFromOldContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateTimeSchemeVariablesFromOldContributions.md)
- [UpdateAuxiliaryVariable](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable.md)
- [UpdateEntityAuxiliaryVariableContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityAuxiliaryVariableContributions.md)
- [CheckVariables](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CheckVariables.md)
- [CalculateBossakConstants](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateBossakConstants.md)
- [GetTimeStep](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GetTimeStep.md)
- [Hash](../../../classes/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/Hash.md)
- [Pred](../../../classes/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/Pred.md)
- [GatherVariables](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GatherVariables.md)
- [SetToZero_AdjointVars](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/SetToZero_AdjointVars.md)
- [Assemble_AdjointVars](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Assemble_AdjointVars.md)

# Imports

- `<vector>`
- `<string>`
- `<unordered_set>`
- `<functional>`
- `"includes/define.h"`
- `"includes/checks.h"`
- `"includes/kratos_parameters.h"`
- `"solving_strategies/schemes/scheme.h"`
- `"response_functions/adjoint_response_function.h"`
- `"utilities/variable_utils.h"`
- `"utilities/indirect_scalar.h"`
- `"utilities/adjoint_extensions.h"`
- `"utilities/atomic_utilities.h"`
- `"utilities/parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)