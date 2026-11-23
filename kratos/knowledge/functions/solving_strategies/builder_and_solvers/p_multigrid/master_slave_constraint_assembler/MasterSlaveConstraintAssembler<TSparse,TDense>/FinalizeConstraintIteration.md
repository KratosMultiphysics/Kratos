---
type: C++ Method
title: FinalizeConstraintIteration
resource: solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler.cpp#L451-L477
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

`bool MasterSlaveConstraintAssembler<TSparse,TDense>::FinalizeConstraintIteration(typename TSparse::MatrixType& rLhs, typename TSparse::VectorType& rSolution, typename TSparse::VectorType& rRhs, typename Base::DofSet::iterator itDofBegin, typename Base::DofSet::iterator itDofEnd, PMGStatusStream::Report& rReport, PMGStatusStream& rStream)`

# Calls

- [Submit](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit.md)
- [Tag](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/Report/Tag.md)