---
type: C++ Method
title: BuildRHSNoDirichletWithoutConstraints
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L1823-L1879
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleRHSWithoutConstraints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSNoDirichlet
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BuildRHSNoDirichletWithoutConstraints( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rb )`

# Calls

- [AssembleRHSWithoutConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleRHSWithoutConstraints.md)

# Called by

- [BuildRHSNoDirichlet](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSNoDirichlet.md)