---
type: C++ Method
title: CalculateSecondDerivativeContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L470-L478
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntitySecondDerivativeContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateSecondDerivativeContributions( Element& rElement, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [CalculateEntitySecondDerivativeContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntitySecondDerivativeContributions.md)