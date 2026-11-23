---
type: C++ Method
title: Assemble
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L1543-L1562
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

`void Assemble( TSystemMatrixType& A, TSystemVectorType& b, const LocalSystemMatrixType& LHS_Contribution, const LocalSystemVectorType& RHS_Contribution, Element::EquationIdVectorType& EquationId )`

# Calls

- [AssembleRowContribution](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/AssembleRowContribution.md)