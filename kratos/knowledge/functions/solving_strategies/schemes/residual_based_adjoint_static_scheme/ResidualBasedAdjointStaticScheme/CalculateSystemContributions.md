---
type: C++ Method
title: CalculateSystemContributions
resource: solving_strategies/schemes/residual_based_adjoint_static_scheme.h#L141-L170
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

`void CalculateSystemContributions(Element& rCurrentElement, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHSContribution, Element::EquationIdVectorType& rEquationId, const ProcessInfo& rCurrentProcessInfo) override`

# Calls

- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)