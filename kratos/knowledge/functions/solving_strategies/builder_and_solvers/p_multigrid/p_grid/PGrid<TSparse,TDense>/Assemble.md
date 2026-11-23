---
type: C++ Method
title: Assemble
resource: solving_strategies/builder_and_solvers/p_multigrid/p_grid.cpp#L165-L297
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler/ConstraintAssembler/GetImposition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PGrid<TSparse,TDense>::Assemble(ModelPart& rModelPart, const typename TParentSparse::MatrixType* pParentLhs, const typename TParentSparse::VectorType* pParentRhs, const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler, IndirectDofSet& rParentDofSet)`

# Calls

- [MatrixMultiplication](../../../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication.md)
- [TransposeMatrix](../../../../../../functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix.md)
- [GetImposition](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler/ConstraintAssembler/GetImposition.md)