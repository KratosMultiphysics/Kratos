---
type: C++ Method
title: Predict
resource: solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp#L225-L319
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
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/UpdateVelocity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/UpdateAcceleration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Predict( ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [HasDofFor](../../../../../functions/includes/node/HasDofFor.md)
- [GetDofPosition](../../../../../functions/includes/node/GetDofPosition.md)
- [UpdateVelocity](../../../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/UpdateVelocity.md)
- [UpdateAcceleration](../../../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/UpdateAcceleration.md)