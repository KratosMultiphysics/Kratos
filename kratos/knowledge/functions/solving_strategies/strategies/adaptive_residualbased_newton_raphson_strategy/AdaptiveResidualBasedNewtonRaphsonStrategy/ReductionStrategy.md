---
type: C++ Method
title: ReductionStrategy
resource: solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h#L1077-L1089
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/ReduceTimeStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReductionStrategy(ModelPart& model_part, double NewTime)`

# Calls

- [ReduceTimeStep](../../../../../functions/sources/model_part/ModelPart/ReduceTimeStep.md)
- [X0](../../../../../functions/includes/node/X0.md)
- [Y0](../../../../../functions/includes/node/Y0.md)
- [Z0](../../../../../functions/includes/node/Z0.md)

# Called by

- [Solve](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Solve.md)