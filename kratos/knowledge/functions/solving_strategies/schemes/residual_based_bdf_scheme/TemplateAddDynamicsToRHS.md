---
type: C++ Function
title: TemplateAddDynamicsToRHS
resource: solving_strategies/schemes/residual_based_bdf_scheme.h#L453-L475
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

`void TemplateAddDynamicsToRHS( TObjectType& rObject, LocalSystemVectorType& rRHS_Contribution, LocalSystemMatrixType& rD, LocalSystemMatrixType& rM, const ProcessInfo& rCurrentProcessInfo )`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)