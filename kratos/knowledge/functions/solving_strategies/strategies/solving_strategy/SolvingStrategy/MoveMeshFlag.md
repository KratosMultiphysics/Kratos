---
type: C++ Method
title: MoveMeshFlag
resource: solving_strategies/strategies/solving_strategy.h#L290-L293
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/CyclicSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/arc_length_strategy/ArcLengthStrategy/SolveSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/SolveSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/SolveSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool MoveMeshFlag()`

# Called by

- [Predict](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Predict.md)
- [CyclicSolve](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/CyclicSolve.md)
- [SolveSolutionStep](../../../../../functions/solving_strategies/strategies/arc_length_strategy/ArcLengthStrategy/SolveSolutionStep.md)
- [Predict](../../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/Predict.md)
- [SolveSolutionStep](../../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/SolveSolutionStep.md)
- [Predict](../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/Predict.md)
- [SolveSolutionStep](../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/SolveSolutionStep.md)