---
type: C++ Method
title: CalculateSensitivity
resource: utilities/geometrical_sensitivity_utility.cpp#L52-L65
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateDeterminantOfJacobianSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateCofactorOfJacobianSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometricalSensitivityUtility::CalculateSensitivity(ShapeParameter Deriv, double& rDetJ_Deriv, ShapeFunctionsGradientType& rDN_DX_Deriv) const`

# Calls

- [CalculateDeterminantOfJacobianSensitivity](../../../../functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateDeterminantOfJacobianSensitivity.md)
- [CalculateCofactorOfJacobianSensitivity](../../../../functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateCofactorOfJacobianSensitivity.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)