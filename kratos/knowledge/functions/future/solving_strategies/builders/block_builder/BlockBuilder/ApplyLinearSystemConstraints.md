---
type: C++ Method
title: ApplyLinearSystemConstraints
resource: future/solving_strategies/builders/block_builder.h#L172-L183
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyLinearSystemConstraints(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override`

# Calls

- [ApplyBlockBuildMasterSlaveConstraints](../../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildMasterSlaveConstraints.md)
- [pGetVector](../../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)