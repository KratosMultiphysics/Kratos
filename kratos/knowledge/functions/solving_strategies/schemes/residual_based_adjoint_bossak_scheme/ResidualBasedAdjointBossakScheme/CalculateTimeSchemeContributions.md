---
type: C++ Method
title: CalculateTimeSchemeContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L542-L553
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateTimeSchemeContributions( Element& rElement, LocalSystemVectorType& rAdjointTimeSchemeValues2, LocalSystemVectorType& rAdjointTimeSchemeValues3, AdjointResponseFunction& rAdjointResponseFunction, const BossakConstants& rBossakConstants, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [CalculateEntityTimeSchemeContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions.md)