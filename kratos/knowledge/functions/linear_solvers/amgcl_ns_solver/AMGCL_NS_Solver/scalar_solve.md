---
type: C++ Method
title: scalar_solve
resource: linear_solvers/amgcl_ns_solver.h#L203-L233
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/amgcl_ns_solver/AMGCL_NS_Solver/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::tuple<size_t, double> scalar_solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) const`

# Called by

- [Solve](../../../../functions/linear_solvers/amgcl_ns_solver/AMGCL_NS_Solver/Solve.md)