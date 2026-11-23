---
type: C++ Method
title: ConstructMatrixStructureWithConstraints
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L777-L942
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
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ConstructMatrixStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void ConstructMatrixStructureWithConstraints( typename TSchemeType::Pointer pScheme, TSystemMatrixType& rA, ModelPart& rModelPart )`

# Calls

- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [ConstructMatrixStructure](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ConstructMatrixStructure.md)