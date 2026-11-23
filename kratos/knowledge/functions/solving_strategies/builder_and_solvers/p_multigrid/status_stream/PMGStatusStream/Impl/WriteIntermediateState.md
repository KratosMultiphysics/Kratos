---
type: C++ Method
title: WriteIntermediateState
resource: solving_strategies/builder_and_solvers/p_multigrid/status_stream.cpp#L197-L253
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/MakeVtuOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ProjectGrid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteIntermediateState(const PMGStatusStream::Report& rReport)`

# Calls

- [visit](../../../../../../../functions/python/variable_indexing_python/visit.md)
- [MakeVtuOutput](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/MakeVtuOutput.md)
- [ProjectGrid](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ProjectGrid.md)
- [emplace_back](../../../../../../../functions/containers/pointer_vector/emplace_back.md)

# Called by

- [Submit](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit.md)