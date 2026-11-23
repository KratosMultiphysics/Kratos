---
type: C++ Method
title: ApplyDirichletConditionsRHS
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L2002-L2028
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSWithConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyDirichletConditionsRHS( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rb )`

# Called by

- [BuildRHSWithConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSWithConstraints.md)