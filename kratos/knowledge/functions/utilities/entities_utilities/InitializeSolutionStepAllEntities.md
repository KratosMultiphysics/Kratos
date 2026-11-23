---
type: C++ Function
title: InitializeSolutionStepAllEntities
resource: utilities/entities_utilities.cpp#L210-L219
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeSolutionStepAllEntities(ModelPart& rModelPart)`

# Called by

- [InitializeSolutionStep](../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeSolutionStep.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/InitializeSolutionStep.md)