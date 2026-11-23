---
type: C++ Method
title: Initialize
resource: solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.cpp#L311-L386
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/GetPenaltyFactor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/InPlaceMatrixAdd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Initialize(typename TSparse::MatrixType& rLhs, typename TSparse::VectorType& rSolution, typename TSparse::VectorType& rRhs, [[maybe_unused]] typename Base::DofSet& rDofs)`

# Calls

- [GetPenaltyFactor](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/GetPenaltyFactor.md)
- [nnz](../../../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [MatrixMultiplication](../../../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication.md)
- [InPlaceMatrixAdd](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/InPlaceMatrixAdd.md)