---
type: C++ Method
title: InitializeSolutionStep
resource: solving_strategies/schemes/residual_based_bdf_custom_scheme.h#L227-L275
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/HasDofFor
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
- [Fix](../../../../../functions/includes/node/Fix.md)