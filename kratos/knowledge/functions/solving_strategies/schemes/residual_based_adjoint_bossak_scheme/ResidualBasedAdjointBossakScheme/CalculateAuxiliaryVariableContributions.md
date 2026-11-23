---
type: C++ Method
title: CalculateAuxiliaryVariableContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L583-L592
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateAuxiliaryVariableContributions( Element& rElement, LocalSystemVectorType& rAdjointAuxiliaryValues, AdjointResponseFunction& rAdjointResponseFunction, const BossakConstants& rBossakConstants, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [CalculateEntityAuxiliaryVariableContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions.md)