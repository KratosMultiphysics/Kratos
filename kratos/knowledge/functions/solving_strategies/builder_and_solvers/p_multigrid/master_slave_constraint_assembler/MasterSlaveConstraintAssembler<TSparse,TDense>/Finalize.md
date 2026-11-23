---
type: C++ Method
title: Finalize
resource: solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler.cpp#L481-L503
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/SetEquationId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MasterSlaveConstraintAssembler<TSparse,TDense>::Finalize([[maybe_unused]] typename TSparse::MatrixType& rLhs, [[maybe_unused]] typename TSparse::VectorType& rSolution, [[maybe_unused]] typename TSparse::VectorType& rRhs, [[maybe_unused]] typename Base::DofSet& rDofSet)`

# Calls

- [SetEquationId](../../../../../../functions/includes/dof/Dof/SetEquationId.md)