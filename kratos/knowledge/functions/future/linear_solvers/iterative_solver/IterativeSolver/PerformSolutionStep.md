---
type: C++ Method
title: PerformSolutionStep
resource: future/linear_solvers/iterative_solver.h#L133-L199
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/containers/linear_system_tags/LinearSystemTags/DenseVectorTagFromString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/containers/linear_system_tags/LinearSystemTags/SparseMatrixTagFromString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/containers/linear_system_tags/LinearSystemTags/DenseMatrixTagFromString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool PerformSolutionStep(LinearSystem<TLinearAlgebra>& rLinearSystem) override`

# Calls

- [DenseVectorTagFromString](../../../../../functions/future/containers/linear_system_tags/LinearSystemTags/DenseVectorTagFromString.md)
- [SparseMatrixTagFromString](../../../../../functions/future/containers/linear_system_tags/LinearSystemTags/SparseMatrixTagFromString.md)
- [pGetVector](../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)
- [DenseMatrixTagFromString](../../../../../functions/future/containers/linear_system_tags/LinearSystemTags/DenseMatrixTagFromString.md)