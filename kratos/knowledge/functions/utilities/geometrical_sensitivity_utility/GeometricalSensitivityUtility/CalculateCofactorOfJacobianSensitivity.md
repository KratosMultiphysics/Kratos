---
type: C++ Method
title: CalculateCofactorOfJacobianSensitivity
resource: utilities/geometrical_sensitivity_utility.cpp#L74-L129
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/CofactorMatrix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateSensitivity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GeometricalSensitivityUtility::MatrixType GeometricalSensitivityUtility::CalculateCofactorOfJacobianSensitivity( ShapeParameter Deriv) const`

# Calls

- [CofactorMatrix](../../../../functions/utilities/math_utils/CofactorMatrix.md)

# Called by

- [CalculateSensitivity](../../../../functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateSensitivity.md)