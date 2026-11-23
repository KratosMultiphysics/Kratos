---
type: C++ Method
title: AssembleRHSWithoutConstraints
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L1913-L1955
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSNoDirichletWithoutConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleWithoutConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssembleRHSWithoutConstraints( TSystemVectorType& rb, const LocalSystemVectorType& rRHSContribution, const Element::EquationIdVectorType& rEquationId )`

# Calls

- [IsFree](../../../../../functions/includes/dof/Dof/IsFree.md)

# Called by

- [BuildRHSNoDirichletWithoutConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSNoDirichletWithoutConstraints.md)
- [AssembleWithoutConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleWithoutConstraints.md)