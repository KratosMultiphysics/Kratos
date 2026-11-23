---
type: C++ Function
title: MatrixMultiplication
resource: utilities/sparse_matrix_multiplication_utility.h#L109-L123
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationSaad
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AllocateSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Initialize
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
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_sparse_matrix_mutiplication_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void MatrixMultiplication( const AMatrix& rA, const BMatrix& rB, CMatrix& rC )`

# Calls

- [MatrixMultiplicationRMerge](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge.md)
- [MatrixMultiplicationSaad](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationSaad.md)

# Called by

- [AllocateSystem](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AllocateSystem.md)
- [Initialize](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Initialize.md)
- [GetIndependentTransform](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Impl/GetIndependentTransform.md)
- [GetDependentTransform](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Impl/GetDependentTransform.md)
- [Initialize](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Initialize.md)
- [Assemble](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Assemble.md)
- [ApplyConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyConstraints.md)
- [BuildWithConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildWithConstraints.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_sparse_matrix_mutiplication_utilities/KRATOS_TEST_CASE_IN_SUITE.md)