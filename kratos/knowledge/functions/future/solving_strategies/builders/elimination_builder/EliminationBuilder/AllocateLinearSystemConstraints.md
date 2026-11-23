---
type: C++ Method
title: AllocateLinearSystemConstraints
resource: future/solving_strategies/builders/elimination_builder.h#L145-L188
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/solving_strategies/builders/builder/Builder/SetUpMasterSlaveConstraintsGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetConstraintsQ
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetConstraintsT
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/EffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AllocateLinearSystemConstraints(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override`

# Calls

- [SetUpMasterSlaveConstraintsGraph](../../../../../../functions/future/solving_strategies/builders/builder/Builder/SetUpMasterSlaveConstraintsGraph.md)
- [pSetConstraintsQ](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetConstraintsQ.md)
- [pSetConstraintsT](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy_data/ImplicitStrategyData/pSetConstraintsT.md)
- [IsFree](../../../../../../functions/includes/dof/Dof/IsFree.md)
- [EffectiveEquationId](../../../../../../functions/includes/dof/Dof/EffectiveEquationId.md)