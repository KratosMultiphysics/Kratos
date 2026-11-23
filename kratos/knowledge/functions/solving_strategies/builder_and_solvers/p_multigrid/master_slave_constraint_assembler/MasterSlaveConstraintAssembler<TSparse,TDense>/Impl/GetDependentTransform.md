---
type: C++ Method
title: GetDependentTransform
resource: solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler.cpp#L80-L101
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/ComputeDependentResidual
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::pair<typename TSparse::MatrixType&, LinearSolver<TSparse,TDense>&> GetDependentTransform()`

# Calls

- [TransposeMatrix](../../../../../../../functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix.md)
- [MatrixMultiplication](../../../../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication.md)

# Called by

- [ComputeDependentResidual](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/ComputeDependentResidual.md)