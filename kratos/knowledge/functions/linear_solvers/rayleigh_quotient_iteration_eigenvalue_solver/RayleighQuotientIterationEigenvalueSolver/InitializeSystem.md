---
type: C++ Method
title: InitializeSystem
resource: linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver.h#L189-L207
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void InitializeSystem( VectorType& rR, const SparseMatrixType& rM )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [Solve](../../../../functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/Solve.md)