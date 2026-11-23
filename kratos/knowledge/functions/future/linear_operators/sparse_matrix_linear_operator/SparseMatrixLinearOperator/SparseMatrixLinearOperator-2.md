---
type: C++ Method
title: SparseMatrixLinearOperator
resource: future/linear_operators/sparse_matrix_linear_operator.h#L77-L81
generated:
  by: okf-rs/0.7.0
---

# Signature

`SparseMatrixLinearOperator(typename MatrixType::Pointer pA) : LinearOperator<TLinearAlgebra>(std::make_pair(pA->size1(), pA->size2())) , mpCsrMatrix(pA)`