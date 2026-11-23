---
type: C++ Method
title: PerformSolutionStep
resource: future/linear_solvers/skyline_lu_factorization_solver.h#L549-L577
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/skyline_lu_factorization_solver/LUSkylineFactorization/copyFromCSRMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool PerformSolutionStep(LinearSystemType& rLinearSystem) override`

# Calls

- [pGetVector](../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)
- [copyFromCSRMatrix](../../../../../functions/future/linear_solvers/skyline_lu_factorization_solver/LUSkylineFactorization/copyFromCSRMatrix.md)