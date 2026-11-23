---
type: C++ Method
title: FinalizeConstraintIteration
resource: solving_strategies/builder_and_solvers/p_multigrid/noop_constraint_assembler.hpp#L45-L57
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

`bool FinalizeConstraintIteration(typename TSparse::MatrixType& rLhs, typename TSparse::VectorType& rSolution, typename TSparse::VectorType& rRhs, typename Base::DofSet::iterator itDofBegin, typename Base::DofSet::iterator itDofEnd, PMGStatusStream::Report& rReport, PMGStatusStream& rStream) override`

# Calls

- [Submit](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit.md)
- [Tag](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/Report/Tag.md)