---
type: C++ Method
title: ApplyBlockBuildDirichletConditions
resource: future/solving_strategies/builders/block_builder.h#L266-L297
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/EffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/GetDiagonalScalingFactor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyBlockBuildDirichletConditions( const DofsArrayType& rDofArray, MatrixType& rLHS, VectorType& rRHS) const`

# Calls

- [EffectiveEquationId](../../../../../../functions/includes/dof/Dof/EffectiveEquationId.md)
- [GetDiagonalScalingFactor](../../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/GetDiagonalScalingFactor.md)