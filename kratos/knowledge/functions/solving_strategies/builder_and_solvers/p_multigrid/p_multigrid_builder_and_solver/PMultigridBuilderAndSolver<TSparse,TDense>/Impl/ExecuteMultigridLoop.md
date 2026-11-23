---
type: C++ Method
title: ExecuteMultigridLoop
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L103-L160
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
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

`void ExecuteMultigridLoop(typename Interface::TSystemMatrixType& rLhs, typename Interface::TSystemVectorType& rSolution, const typename Interface::TSystemVectorType& rRhs, typename Interface::TSystemVectorType& rSolutionUpdate, typename Interface::TSystemVectorType& rResidual, const typename TSparse::DataType InitialResidualNorm, PMGStatusStream& rStream, PMGStatusStream::Report& rReport)`

# Calls

- [visit](../../../../../../../functions/python/variable_indexing_python/visit.md)
- [GetLinearSystemSolver](../../../../../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetLinearSystemSolver.md)
- [Submit](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit.md)
- [Tag](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/Report/Tag.md)