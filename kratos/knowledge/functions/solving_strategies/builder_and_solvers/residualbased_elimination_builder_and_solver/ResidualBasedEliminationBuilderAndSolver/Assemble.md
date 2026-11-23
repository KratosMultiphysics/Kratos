---
type: C++ Method
title: Assemble
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h#L949-L982
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/AssembleRowContributionFreeDofs
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Assemble( TSystemMatrixType& rA, TSystemVectorType& rb, const LocalSystemMatrixType& rLHSContribution, const LocalSystemVectorType& rRHSContribution, const Element::EquationIdVectorType& rEquationId #ifdef USE_LOCKS_IN_ASSEMBLY ,std::vector< omp_lock_t >& rLockArray #endif )`

# Calls

- [AssembleRowContributionFreeDofs](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/AssembleRowContributionFreeDofs.md)