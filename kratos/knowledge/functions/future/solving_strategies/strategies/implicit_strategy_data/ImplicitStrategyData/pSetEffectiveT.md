---
type: C++ Method
title: pSetEffectiveT
resource: future/solving_strategies/strategies/implicit_strategy_data.h#L193-L196
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/ApplyLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void pSetEffectiveT(typename MatrixType::Pointer pEffectiveT)`

# Called by

- [ApplyBlockBuildMasterSlaveConstraints](../../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildMasterSlaveConstraints.md)
- [ApplyLinearSystemConstraints](../../../../../../functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/ApplyLinearSystemConstraints.md)