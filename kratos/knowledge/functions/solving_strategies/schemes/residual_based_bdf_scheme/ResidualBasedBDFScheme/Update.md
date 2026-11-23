---
type: C++ Method
title: Update
resource: solving_strategies/schemes/residual_based_bdf_scheme.h#L186-L202
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/dof_updater/DofUpdater/UpdateDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/UpdateDerivatives
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Update( ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [UpdateDofs](../../../../../functions/utilities/dof_updater/DofUpdater/UpdateDofs.md)
- [UpdateDerivatives](../../../../../functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/UpdateDerivatives.md)