---
type: C++ Method
title: Solve
resource: linear_solvers/power_iteration_highest_eigenvalue_solver.h#L178-L246
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/random_initializer_utility/RandomInitializeUtility/RandomInitialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Solve( SparseMatrixType& K, SparseMatrixType& M, DenseVectorType& Eigenvalues, DenseMatrixType& Eigenvectors ) override`

# Calls

- [RandomInitialize](../../../../functions/utilities/random_initializer_utility/RandomInitializeUtility/RandomInitialize.md)