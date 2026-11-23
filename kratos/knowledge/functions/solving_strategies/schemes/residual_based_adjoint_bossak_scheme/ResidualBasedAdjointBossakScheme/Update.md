---
type: C++ Method
title: Update
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L226-L242
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/dof_updater/DofUpdater/UpdateDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Update( ModelPart& rModelPart, DofsArrayType& rDofSet, SystemMatrixType& rA, SystemVectorType& rDx, SystemVectorType& rb) override`

# Calls

- [UpdateDofs](../../../../../functions/utilities/dof_updater/DofUpdater/UpdateDofs.md)
- [UpdateTimeSchemeAdjoints](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints.md)