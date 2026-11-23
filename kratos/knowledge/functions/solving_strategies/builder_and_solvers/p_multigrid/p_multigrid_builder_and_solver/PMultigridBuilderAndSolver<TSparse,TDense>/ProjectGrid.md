---
type: C++ Method
title: ProjectGrid
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L1193-L1283
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid/GetChild
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PMultigridBuilderAndSolver<TSparse,TDense>::ProjectGrid(int GridLevel, const typename TSparse::MatrixType& rRootLhs, const typename TSparse::VectorType& rRootSolution, const typename TSparse::VectorType& rRootRhs)`

# Calls

- [visit](../../../../../../functions/python/variable_indexing_python/visit.md)
- [GetChild](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid/GetChild.md)
- [fill](../../../../../../functions/includes/amatrix_interface/Matrix/fill.md)

# Called by

- [WriteIntermediateState](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState.md)