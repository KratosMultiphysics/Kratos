---
type: C++ Method
title: InitializeConditions
resource: solving_strategies/schemes/scheme.h#L246-L257
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/scheme/Scheme/SetConditionsAreInitialized
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void InitializeConditions(ModelPart& rModelPart)`

# Calls

- [SetConditionsAreInitialized](../../../../../functions/solving_strategies/schemes/scheme/Scheme/SetConditionsAreInitialized.md)

# Called by

- [Initialize](../../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/Initialize.md)
- [Initialize](../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/Initialize.md)