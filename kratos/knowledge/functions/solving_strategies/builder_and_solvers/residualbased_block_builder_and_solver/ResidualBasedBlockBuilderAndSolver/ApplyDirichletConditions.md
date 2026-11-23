---
type: C++ Method
title: ApplyDirichletConditions
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L859-L908
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spaces/ublas_space/CheckAndCorrectZeroDiagonalValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyDirichletConditions( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [CheckAndCorrectZeroDiagonalValues](../../../../../functions/spaces/ublas_space/CheckAndCorrectZeroDiagonalValues.md)