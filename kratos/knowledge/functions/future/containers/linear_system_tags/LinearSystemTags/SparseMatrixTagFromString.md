---
type: C++ Method
title: SparseMatrixTagFromString
resource: future/containers/linear_system_tags.h#L83-L99
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/linear_solvers/iterative_solver/IterativeSolver/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/iterative_solver/IterativeSolver/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/iterative_solver/IterativeSolver/PerformSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/iterative_solver/IterativeSolver/FinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static SparseMatrixTag SparseMatrixTagFromString(const std::string &rTagName)`

# Called by

- [Initialize](../../../../../functions/future/linear_solvers/iterative_solver/IterativeSolver/Initialize.md)
- [InitializeSolutionStep](../../../../../functions/future/linear_solvers/iterative_solver/IterativeSolver/InitializeSolutionStep.md)
- [PerformSolutionStep](../../../../../functions/future/linear_solvers/iterative_solver/IterativeSolver/PerformSolutionStep.md)
- [FinalizeSolutionStep](../../../../../functions/future/linear_solvers/iterative_solver/IterativeSolver/FinalizeSolutionStep.md)