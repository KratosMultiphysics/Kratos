---
type: C++ Method
title: SturmSequenceCheck
resource: linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver.h#L213-L231
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/linear_solvers/skyline_lu_factorization_solver/LUSkylineFactorization/copyFromCSRMatrix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType SturmSequenceCheck(SparseMatrixType& ShiftedK)`

# Calls

- [copyFromCSRMatrix](../../../../functions/future/linear_solvers/skyline_lu_factorization_solver/LUSkylineFactorization/copyFromCSRMatrix.md)

# Called by

- [Solve](../../../../functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/Solve.md)