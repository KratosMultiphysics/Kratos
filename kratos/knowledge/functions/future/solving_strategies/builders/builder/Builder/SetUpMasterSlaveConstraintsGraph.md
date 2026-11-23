---
type: C++ Method
title: SetUpMasterSlaveConstraintsGraph
resource: future/solving_strategies/builders/builder.h#L203-L253
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/EffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/AllocateLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void SetUpMasterSlaveConstraintsGraph( const DofsArrayType& rDofSet, const DofsArrayType& rEffectiveDofSet, SparseGraphType& rConstraintsSparseGraph)`

# Calls

- [EffectiveEquationId](../../../../../../functions/includes/dof/Dof/EffectiveEquationId.md)

# Called by

- [AllocateLinearSystemConstraints](../../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/AllocateLinearSystemConstraints.md)
- [AllocateLinearSystemConstraints](../../../../../../functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystemConstraints.md)