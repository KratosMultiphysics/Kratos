---
type: C++ Method
title: Solve
resource: linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver.h#L240-L390
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/InitializeSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/skyline_lu_factorization_solver/LUSkylineFactorization/copyFromCSRMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/SturmSequenceCheck
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Solve( SparseMatrixType& K, SparseMatrixType& M, DenseVectorType& Eigenvalues, DenseMatrixType& Eigenvectors ) override`

# Calls

- [InitializeSystem](../../../../functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/InitializeSystem.md)
- [copyFromCSRMatrix](../../../../functions/future/linear_solvers/skyline_lu_factorization_solver/LUSkylineFactorization/copyFromCSRMatrix.md)
- [SturmSequenceCheck](../../../../functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/SturmSequenceCheck.md)