---
type: C++ Method
title: Initialize
resource: solving_strategies/strategies/explicit_solving_strategy.h#L171-L189
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/entities_utilities/InitializeAllEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetResetDofSetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/GetMoveMeshFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetResetLumpedMassVectorFlag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Initialize() override`

# Calls

- [InitializeAllEntities](../../../../../functions/utilities/entities_utilities/InitializeAllEntities.md)
- [SetResetDofSetFlag](../../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetResetDofSetFlag.md)
- [GetMoveMeshFlag](../../../../../functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/GetMoveMeshFlag.md)
- [SetResetLumpedMassVectorFlag](../../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetResetLumpedMassVectorFlag.md)