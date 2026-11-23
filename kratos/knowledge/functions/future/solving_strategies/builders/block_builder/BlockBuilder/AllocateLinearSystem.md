---
type: C++ Method
title: AllocateLinearSystem
resource: future/solving_strategies/builders/block_builder.h#L120-L148
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetLinearSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/RequiresEffectiveDofSet
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
- [RequiresEffectiveDofSet](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/RequiresEffectiveDofSet.md)
- [pSetEffectiveLinearSystem](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetEffectiveLinearSystem.md)