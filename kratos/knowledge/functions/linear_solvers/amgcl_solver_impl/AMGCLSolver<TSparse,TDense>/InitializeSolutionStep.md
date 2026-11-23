---
type: C++ Method
title: InitializeSolutionStep
resource: linear_solvers/amgcl_solver_impl.hpp#L454-L679
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AMGCLSolver<TSparse,TDense>::InitializeSolutionStep(SparseMatrixType& rLhs, VectorType& rSolution, VectorType& rRhs)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)