---
type: C++ Method
title: pSetConstraintsQ
resource: future/solving_strategies/strategies/implicit_strategy_data.h#L247-L250
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/AllocateLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void pSetConstraintsQ(typename VectorType::Pointer pConstraintsQ)`

# Called by

- [AllocateLinearSystemConstraints](../../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/AllocateLinearSystemConstraints.md)
- [AllocateLinearSystemConstraints](../../../../../../functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystemConstraints.md)