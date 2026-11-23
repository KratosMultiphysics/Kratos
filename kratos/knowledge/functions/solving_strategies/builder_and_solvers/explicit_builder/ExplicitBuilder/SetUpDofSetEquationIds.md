---
type: C++ Method
title: SetUpDofSetEquationIds
resource: solving_strategies/builder_and_solvers/explicit_builder.h#L664-L677
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/SetEquationId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void SetUpDofSetEquationIds()`

# Calls

- [SetEquationId](../../../../../functions/includes/dof/Dof/SetEquationId.md)

# Called by

- [Initialize](../../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/Initialize.md)
- [InitializeSolutionStep](../../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/InitializeSolutionStep.md)