---
type: C++ Method
title: CalculateGradientContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L430-L438
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityGradientContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateGradientContributions( Element& rElement, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [CalculateEntityGradientContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityGradientContributions.md)