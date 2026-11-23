---
type: C++ Method
title: ExecuteConstraintLoop
resource: solving_strategies/builder_and_solvers/p_multigrid/p_grid.cpp#L409-L455
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/Report/Tag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PGrid<TSparse,TDense>::ExecuteConstraintLoop(PMGStatusStream& rStream, PMGStatusStream::Report& rReport)`

# Calls

- [Submit](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit.md)
- [Tag](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/Report/Tag.md)