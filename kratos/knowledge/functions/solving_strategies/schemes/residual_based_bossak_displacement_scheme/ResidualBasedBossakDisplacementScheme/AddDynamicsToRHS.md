---
type: C++ Method
title: AddDynamicsToRHS
resource: solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp#L614-L642
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

`void AddDynamicsToRHS( Element& rElement, LocalSystemVectorType& RHS_Contribution, LocalSystemMatrixType& D, LocalSystemMatrixType& M, const ProcessInfo& rCurrentProcessInfo ) override`

# Calls

- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)