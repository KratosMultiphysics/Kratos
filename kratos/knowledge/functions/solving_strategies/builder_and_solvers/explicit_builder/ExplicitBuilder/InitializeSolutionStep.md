---
type: C++ Method
title: InitializeSolutionStep
resource: solving_strategies/builder_and_solvers/explicit_builder.h#L387-L403
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpDofSetEquationIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpLumpedMassVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/InitializeDofSetReactions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void InitializeSolutionStep(ModelPart& rModelPart)`

# Calls

- [SetUpDofSetEquationIds](../../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpDofSetEquationIds.md)
- [SetUpLumpedMassVector](../../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpLumpedMassVector.md)
- [InitializeDofSetReactions](../../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/InitializeDofSetReactions.md)