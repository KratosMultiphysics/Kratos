---
type: C++ Method
title: ApplyDirichletConditions
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier.h#L402-L484
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spaces/ublas_space/GetScaleNorm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyDirichletConditions( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [GetScaleNorm](../../../../../functions/spaces/ublas_space/GetScaleNorm.md)