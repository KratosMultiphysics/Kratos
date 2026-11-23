---
type: C++ Method
title: SolveSolutionStep
resource: solving_strategies/strategies/explicit_solving_strategy.h#L230-L252
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/entities_utilities/InitializeNonLinearIterationAllEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/GetMoveMeshFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/entities_utilities/FinalizeNonLinearIterationAllEntities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SolveSolutionStep() override`

# Calls

- [InitializeNonLinearIterationAllEntities](../../../../../functions/utilities/entities_utilities/InitializeNonLinearIterationAllEntities.md)
- [GetMoveMeshFlag](../../../../../functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/GetMoveMeshFlag.md)
- [FinalizeNonLinearIterationAllEntities](../../../../../functions/utilities/entities_utilities/FinalizeNonLinearIterationAllEntities.md)