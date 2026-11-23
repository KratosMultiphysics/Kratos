---
type: C++ Method
title: CalculateRHSContribution
resource: solving_strategies/schemes/residual_based_adjoint_static_scheme.h#L185-L211
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateRHSContribution( Element& rCurrentElement, LocalSystemVectorType& rRHSContribution, Element::EquationIdVectorType& rEquationId, const ProcessInfo& rCurrentProcessInfo) override`

# Calls

- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)