---
type: C++ Method
title: ApplyConstraints
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L954-L995
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/GetScaleNorm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyConstraints( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rb ) override`

# Calls

- [MatrixMultiplication](../../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication.md)
- [GetScaleNorm](../../../../../functions/spaces/ublas_space/GetScaleNorm.md)