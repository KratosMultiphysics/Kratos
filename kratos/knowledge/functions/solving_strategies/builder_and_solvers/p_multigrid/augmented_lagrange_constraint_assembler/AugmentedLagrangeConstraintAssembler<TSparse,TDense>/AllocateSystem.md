---
type: C++ Method
title: AllocateSystem
resource: solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.cpp#L241-L266
visibility: private
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
---

# Signature

`void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AllocateSystem(typename TSparse::MatrixType& rLhs, typename TSparse::VectorType& rSolution, typename TSparse::VectorType& rRhs, typename Base::DofSet& rDofs)`

# Calls

- [TransposeMatrix](../../../../../../functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix.md)
- [MatrixMultiplication](../../../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication.md)