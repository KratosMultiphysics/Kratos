---
type: C++ Method
title: UpdateConstraintsOnlyDofs
resource: future/solving_strategies/schemes/implicit_scheme.h#L1119-L1147
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/EffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Update
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void UpdateConstraintsOnlyDofs( const VectorType& rEffectiveDx, DofsArrayType& rDofSet, DofsArrayType& rEffectiveDofSet)`

# Calls

- [IsFree](../../../../../../functions/includes/dof/Dof/IsFree.md)
- [EffectiveEquationId](../../../../../../functions/includes/dof/Dof/EffectiveEquationId.md)

# Called by

- [Update](../../../../../../functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Update.md)