---
type: C++ Method
title: SetUpSystemWithConstraints
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L1497-L1606
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/SetEquationId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetUpSystemWithConstraints(ModelPart& rModelPart)`

# Calls

- [SetEquationId](../../../../../functions/includes/dof/Dof/SetEquationId.md)

# Called by

- [SetUpSystem](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpSystem.md)