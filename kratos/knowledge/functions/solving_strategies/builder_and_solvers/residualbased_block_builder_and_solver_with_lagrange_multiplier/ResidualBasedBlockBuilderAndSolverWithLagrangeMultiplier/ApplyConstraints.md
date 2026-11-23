---
type: C++ Method
title: ApplyConstraints
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier.h#L564-L701
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spaces/ublas_space/GetDiagonalNorm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/GetAveragevalueDiagonal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/AssembleSparseMatrixByBlocks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ComputeRHSLMContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyConstraints( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rb ) override`

# Calls

- [GetDiagonalNorm](../../../../../functions/spaces/ublas_space/GetDiagonalNorm.md)
- [GetAveragevalueDiagonal](../../../../../functions/spaces/ublas_space/GetAveragevalueDiagonal.md)
- [nnz](../../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [AssembleSparseMatrixByBlocks](../../../../../functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/AssembleSparseMatrixByBlocks.md)
- [ComputeRHSLMContributions](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ComputeRHSLMContributions.md)