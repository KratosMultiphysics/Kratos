---
type: C++ Method
title: pSetEffectiveLinearSystem
resource: future/solving_strategies/strategies/implicit_strategy_data.h#L166-L169
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/AllocateLinearSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void pSetEffectiveLinearSystem(typename LinearSystemType::Pointer pEffectiveLinearSystem)`

# Called by

- [AllocateLinearSystem](../../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/AllocateLinearSystem.md)
- [AllocateLinearSystem](../../../../../../functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystem.md)