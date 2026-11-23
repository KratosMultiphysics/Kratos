---
type: C++ Method
title: CalculateSolutionVector
resource: future/solving_strategies/builders/builder.h#L281-L312
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/EffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateSolutionVector( const DofsArrayType& rEffectiveDofSet, const MatrixType& rConstraintsRelationMatrix, const VectorType& rConstraintsConstantVector, VectorType& rSolutionVector) const`

# Calls

- [EffectiveEquationId](../../../../../../functions/includes/dof/Dof/EffectiveEquationId.md)

# Called by

- [Predict](../../../../../../functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Predict.md)