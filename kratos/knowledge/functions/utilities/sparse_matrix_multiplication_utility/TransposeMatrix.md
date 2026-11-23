---
type: C++ Function
title: TransposeMatrix
resource: utilities/sparse_matrix_multiplication_utility.h#L489-L561
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sparse_matrix_multiplication_utility/SortRows
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/CreateSolutionMatrix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AllocateSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Impl/GetIndependentTransform
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Impl/GetDependentTransform
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Assemble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ApplyDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixTranspose/__sparse_matrix_transpose
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void TransposeMatrix( AMatrix& rA, const BMatrix& rB, const double Factor = 1.0 )`

# Calls

- [SortRows](../../../functions/utilities/sparse_matrix_multiplication_utility/SortRows.md)
- [CreateSolutionMatrix](../../../functions/utilities/sparse_matrix_multiplication_utility/CreateSolutionMatrix.md)

# Called by

- [AllocateSystem](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AllocateSystem.md)
- [GetIndependentTransform](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Impl/GetIndependentTransform.md)
- [GetDependentTransform](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Impl/GetDependentTransform.md)
- [Initialize](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Initialize.md)
- [Assemble](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Assemble.md)
- [ApplyDirichletConditions](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ApplyDirichletConditions.md)
- [__sparse_matrix_transpose](../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixTranspose/__sparse_matrix_transpose.md)