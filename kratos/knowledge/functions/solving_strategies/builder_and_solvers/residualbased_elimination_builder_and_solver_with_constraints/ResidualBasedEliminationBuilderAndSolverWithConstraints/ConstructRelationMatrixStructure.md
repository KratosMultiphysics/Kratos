---
type: C++ Method
title: ConstructRelationMatrixStructure
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L950-L1055
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ResizeAndInitializeVectors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void ConstructRelationMatrixStructure( typename TSchemeType::Pointer pScheme, TSystemMatrixType& rT, ModelPart& rModelPart )`

# Calls

- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [ResizeAndInitializeVectors](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ResizeAndInitializeVectors.md)