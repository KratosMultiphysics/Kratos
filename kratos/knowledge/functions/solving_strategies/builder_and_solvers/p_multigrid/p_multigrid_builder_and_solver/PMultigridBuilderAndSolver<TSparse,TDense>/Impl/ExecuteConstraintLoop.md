---
type: C++ Method
title: ExecuteConstraintLoop
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L163-L244
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/builder_and_solver/GetLinearSystemSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/Report/Tag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`[[nodiscard]] PMGStatusStream::Report ExecuteConstraintLoop(typename Interface::TSystemMatrixType& rLhs, typename Interface::TSystemVectorType& rSolution, typename Interface::TSystemVectorType& rRhs, PMGStatusStream& rStream)`

# Calls

- [GetLinearSystemSolver](../../../../../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetLinearSystemSolver.md)
- [Submit](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit.md)
- [Tag](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/Report/Tag.md)