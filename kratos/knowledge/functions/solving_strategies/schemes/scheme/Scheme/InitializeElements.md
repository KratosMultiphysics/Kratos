---
type: C++ Method
title: InitializeElements
resource: solving_strategies/schemes/scheme.h#L230-L239
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/scheme/Scheme/SetElementsAreInitialized
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void InitializeElements( ModelPart& rModelPart)`

# Calls

- [SetElementsAreInitialized](../../../../../functions/solving_strategies/schemes/scheme/Scheme/SetElementsAreInitialized.md)

# Called by

- [Initialize](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Initialize.md)
- [Initialize](../../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/Initialize.md)
- [Initialize](../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/Initialize.md)