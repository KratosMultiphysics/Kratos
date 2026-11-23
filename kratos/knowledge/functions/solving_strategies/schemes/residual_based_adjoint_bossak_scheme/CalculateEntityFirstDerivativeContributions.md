---
type: C++ Function
title: CalculateEntityFirstDerivativeContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L735-L749
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateFirstDerivativeContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEntityFirstDerivativeContributions( TEntityType& rCurrentEntity, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)

# Called by

- [CalculateFirstDerivativeContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateFirstDerivativeContributions.md)