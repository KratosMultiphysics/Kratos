---
type: C++ Method
title: CalculateFirstDerivativeContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L450-L458
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityFirstDerivativeContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateFirstDerivativeContributions( Element& rElement, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [CalculateEntityFirstDerivativeContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityFirstDerivativeContributions.md)