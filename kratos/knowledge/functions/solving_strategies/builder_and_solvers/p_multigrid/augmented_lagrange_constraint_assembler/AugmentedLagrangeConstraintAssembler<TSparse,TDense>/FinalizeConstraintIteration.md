---
type: C++ Method
title: FinalizeConstraintIteration
resource: solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.cpp#L417-L456
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

`bool AugmentedLagrangeConstraintAssembler<TSparse,TDense>::FinalizeConstraintIteration(typename TSparse::MatrixType& rLhs, typename TSparse::VectorType& rSolution, [[maybe_unused]] typename TSparse::VectorType& rRhs, [[maybe_unused]] typename Base::DofSet::iterator itDofBegin, [[maybe_unused]] typename Base::DofSet::iterator itDofEnd, PMGStatusStream::Report& rReport, PMGStatusStream& rStream)`

# Calls

- [Submit](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Submit.md)
- [Tag](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/Report/Tag.md)