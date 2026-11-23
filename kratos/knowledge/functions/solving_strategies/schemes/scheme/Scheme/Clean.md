---
type: C++ Method
title: Clean
resource: solving_strategies/schemes/scheme.h#L489-L493
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Solve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/FinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/FinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void Clean()`

# Called by

- [Solve](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Solve.md)
- [FinalizeSolutionStep](../../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/FinalizeSolutionStep.md)
- [FinalizeSolutionStep](../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/FinalizeSolutionStep.md)