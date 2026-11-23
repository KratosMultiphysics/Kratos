---
type: C++ Method
title: AllocateLinearSystem
resource: future/solving_strategies/builders/elimination_builder.h#L114-L143
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetLinearSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetEffectiveLinearSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AllocateLinearSystem( const SparseGraphType& rSparseGraph, ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override`

# Calls

- [pSetLinearSystem](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetLinearSystem.md)
- [pSetEffectiveLinearSystem](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetEffectiveLinearSystem.md)