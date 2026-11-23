---
type: C++ Method
title: AddDynamicsToRHS
resource: solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme.h#L404-L423
visibility: private
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

`void AddDynamicsToRHS( Element& rElement, LocalSystemVectorType& rRHSContribution, LocalSystemMatrixType& rD, LocalSystemMatrixType& rM, const ProcessInfo& rCurrentProcessInfo ) override`

# Calls

- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)