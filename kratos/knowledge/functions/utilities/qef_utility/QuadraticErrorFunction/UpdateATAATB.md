---
type: C++ Method
title: UpdateATAATB
resource: utilities/qef_utility.cpp#L113-L121
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/QuadraticErrorFunctionPoint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void QuadraticErrorFunction::UpdateATAATB(AuxiliaryClasses& rAuxiliaryClasses)`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [QuadraticErrorFunctionPoint](../../../../functions/utilities/qef_utility/QuadraticErrorFunction/QuadraticErrorFunctionPoint.md)