---
type: C++ Method
title: CalculateSystemContributions
resource: solving_strategies/schemes/residual_based_implicit_time_scheme.h#L158-L171
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_implicit_time_scheme/TCalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateSystemContributions( Element& rCurrentElement, LocalSystemMatrixType& rLHSContribution, LocalSystemVectorType& rRHSContribution, Element::EquationIdVectorType& rEquationId, const ProcessInfo& rCurrentProcessInfo ) override`

# Calls

- [TCalculateSystemContributions](../../../../../functions/solving_strategies/schemes/residual_based_implicit_time_scheme/TCalculateSystemContributions.md)