---
type: C++ Method
title: UpdateVelocity
resource: solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp#L557-L565
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void UpdateVelocity( array_1d<double, 3>& rCurrentVelocity, const array_1d<double, 3>& rDeltaDisplacement, const array_1d<double, 3>& rPreviousVelocity, const array_1d<double, 3>& rPreviousAcceleration )`

# Called by

- [Update](../../../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Update.md)
- [Predict](../../../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Predict.md)