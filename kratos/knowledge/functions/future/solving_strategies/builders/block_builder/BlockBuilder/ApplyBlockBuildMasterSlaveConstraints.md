---
type: C++ Method
title: ApplyBlockBuildMasterSlaveConstraints
resource: future/solving_strategies/builders/block_builder.h#L201-L257
visibility: private
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
  - target: functions/future/containers/linear_system/LinearSystem/pSetVector
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyBlockBuildMasterSlaveConstraints(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)`

# Calls

- [pGetVector](../../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)
- [pSetEffectiveT](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetEffectiveT.md)
- [pSetVector](../../../../../../functions/future/containers/linear_system/LinearSystem/pSetVector.md)

# Called by

- [ApplyLinearSystemConstraints](../../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyLinearSystemConstraints.md)