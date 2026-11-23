---
type: C++ Method
title: ApplyLinearSystemConstraints
resource: future/solving_strategies/builders/elimination_builder.h#L191-L241
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetEffectiveT
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyLinearSystemConstraints(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override`

# Calls

- [pGetVector](../../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)
- [pSetEffectiveT](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetEffectiveT.md)