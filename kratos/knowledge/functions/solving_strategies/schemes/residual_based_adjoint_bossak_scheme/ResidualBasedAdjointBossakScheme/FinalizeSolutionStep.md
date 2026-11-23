---
type: C++ Method
title: FinalizeSolutionStep
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L212-L224
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FinalizeSolutionStep( ModelPart& rModelPart, SystemMatrixType& rA, SystemVectorType& rDx, SystemVectorType& rb) override`

# Calls

- [UpdateAuxiliaryVariable](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable.md)