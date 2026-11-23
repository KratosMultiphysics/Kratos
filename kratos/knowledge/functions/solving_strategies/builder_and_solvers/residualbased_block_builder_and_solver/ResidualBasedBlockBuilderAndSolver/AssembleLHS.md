---
type: C++ Method
title: AssembleLHS
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L1566-L1578
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/AssembleRowContribution
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssembleLHS( TSystemMatrixType& rA, const LocalSystemMatrixType& rLHSContribution, Element::EquationIdVectorType& rEquationId )`

# Calls

- [AssembleRowContribution](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/AssembleRowContribution.md)