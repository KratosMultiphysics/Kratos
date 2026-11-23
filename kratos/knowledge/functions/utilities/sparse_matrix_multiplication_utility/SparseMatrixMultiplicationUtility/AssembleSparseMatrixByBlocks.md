---
type: C++ Method
title: AssembleSparseMatrixByBlocks
resource: utilities/sparse_matrix_multiplication_utility.h#L655-L838
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/ComputeNonZeroBlocks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/check_non_zero_blocks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/ComputeAuxiliarValuesBlocks
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_sparse_matrix_mutiplication_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void AssembleSparseMatrixByBlocks( CompressedMatrix& rMatrix, const DenseMatrix<CompressedMatrix*>& rMatricespBlocks, DenseMatrix<double> ContributionCoefficients = DenseMatrix<double>(0,0), DenseMatrix<bool> TransposeBlocks = DenseMatrix<bool>(0,0) )`

# Calls

- [nnz](../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [ComputeNonZeroBlocks](../../../../functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/ComputeNonZeroBlocks.md)
- [check_non_zero_blocks](../../../../functions/utilities/sparse_matrix_multiplication_utility/check_non_zero_blocks.md)
- [ComputeAuxiliarValuesBlocks](../../../../functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/ComputeAuxiliarValuesBlocks.md)

# Called by

- [ApplyConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_sparse_matrix_mutiplication_utilities/KRATOS_TEST_CASE_IN_SUITE.md)