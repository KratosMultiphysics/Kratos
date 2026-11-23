---
type: C++ Method
title: ReduceTimeStep
resource: sources/model_part.cpp#L273-L289
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/process_info/ProcessInfo/SetCurrentTime
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/ReductionStrategy
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::ReduceTimeStep(ModelPart& rModelPart, double NewTime)`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [SetCurrentTime](../../../../functions/sources/process_info/ProcessInfo/SetCurrentTime.md)

# Called by

- [ReductionStrategy](../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/ReductionStrategy.md)