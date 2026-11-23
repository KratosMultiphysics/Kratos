---
type: C++ Method
title: Update
resource: solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme.h#L166-L191
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/dof_updater/DofUpdater/UpdateDofs
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Update( ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [UpdateDofs](../../../../../functions/utilities/dof_updater/DofUpdater/UpdateDofs.md)