---
type: C++ Method
title: Predict
resource: solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme.h#L202-L293
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/GetDofPosition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Predict( ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [HasDofFor](../../../../../functions/includes/node/HasDofFor.md)
- [GetDofPosition](../../../../../functions/includes/node/GetDofPosition.md)