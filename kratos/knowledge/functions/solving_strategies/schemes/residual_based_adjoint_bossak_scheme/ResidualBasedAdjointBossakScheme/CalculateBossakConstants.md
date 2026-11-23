---
type: C++ Method
title: CalculateBossakConstants
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L1207-L1222
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static BossakConstants CalculateBossakConstants(double Alpha, double DeltaTime)`

# Called by

- [Predict](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Predict.md)
- [InitializeSolutionStep](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/InitializeSolutionStep.md)