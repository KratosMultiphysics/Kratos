---
type: C++ Function
title: CalculateEntityGradientContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L706-L718
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateGradientContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEntityGradientContributions( TEntityType& rCurrentEntity, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)

# Called by

- [CalculateGradientContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateGradientContributions.md)