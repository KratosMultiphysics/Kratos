---
type: C++ Method
title: SetUpDofSetWithConstraints
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L577-L716
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/HasReaction
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpDofSet
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetUpDofSetWithConstraints( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart )`

# Calls

- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)
- [HasReaction](../../../../../functions/includes/dof/Dof/HasReaction.md)

# Called by

- [SetUpDofSet](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpDofSet.md)