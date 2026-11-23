---
type: C++ Method
title: Update
resource: solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp#L187-L214
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/dof_updater/DofUpdater/UpdateDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/UpdateVelocity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/UpdateAcceleration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Update( ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [UpdateDofs](../../../../../functions/utilities/dof_updater/DofUpdater/UpdateDofs.md)
- [UpdateVelocity](../../../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/UpdateVelocity.md)
- [UpdateAcceleration](../../../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/UpdateAcceleration.md)