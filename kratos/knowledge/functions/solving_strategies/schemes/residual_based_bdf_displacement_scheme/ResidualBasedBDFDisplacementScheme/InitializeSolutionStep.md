---
type: C++ Method
title: InitializeSolutionStep
resource: solving_strategies/schemes/residual_based_bdf_displacement_scheme.h#L164-L217
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
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeSolutionStep( ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [HasDofFor](../../../../../functions/includes/node/HasDofFor.md)
- [GetDofPosition](../../../../../functions/includes/node/GetDofPosition.md)
- [Fix](../../../../../functions/includes/node/Fix.md)