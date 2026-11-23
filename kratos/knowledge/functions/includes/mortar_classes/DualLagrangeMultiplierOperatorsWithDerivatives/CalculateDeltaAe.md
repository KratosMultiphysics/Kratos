---
type: C++ Method
title: CalculateDeltaAe
resource: includes/mortar_classes.h#L1732-L1764
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_frobenius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/CheckConditionNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool CalculateDeltaAe(DerivativeDataType& rDerivativeData)`

# Calls

- [norm_frobenius](../../../../functions/includes/amatrix_interface/norm_frobenius.md)
- [CheckConditionNumber](../../../../functions/utilities/math_utils/CheckConditionNumber.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)