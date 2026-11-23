---
type: C++ Method
title: InitializeSolutionStep
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L194-L210
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateBossakConstants
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GetTimeStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateNodeNeighbourCount
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeSolutionStep( ModelPart& rModelPart, SystemMatrixType& rA, SystemVectorType& rDx, SystemVectorType& rb) override`

# Calls

- [CalculateBossakConstants](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateBossakConstants.md)
- [GetTimeStep](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GetTimeStep.md)
- [CalculateNodeNeighbourCount](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateNodeNeighbourCount.md)