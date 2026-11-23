---
type: C++ Method
title: CalculateResidualLocalContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L502-L510
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityResidualLocalContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateResidualLocalContributions( Element& rCurrentElement, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [CalculateEntityResidualLocalContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityResidualLocalContributions.md)