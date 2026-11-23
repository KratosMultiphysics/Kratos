---
type: C++ Method
title: PerformSolutionStep
resource: linear_solvers/amgcl_solver_impl.hpp#L711-L786
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

`bool AMGCLSolver<TSparse,TDense>::PerformSolutionStep(SparseMatrixType& rLhs, VectorType& rSolution, VectorType& rRhs)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)