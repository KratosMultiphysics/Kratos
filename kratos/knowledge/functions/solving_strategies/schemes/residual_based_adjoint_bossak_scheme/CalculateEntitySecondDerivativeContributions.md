---
type: C++ Function
title: CalculateEntitySecondDerivativeContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L766-L781
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSecondDerivativeContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEntitySecondDerivativeContributions( TEntityType& rCurrentEntity, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)

# Called by

- [CalculateSecondDerivativeContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSecondDerivativeContributions.md)