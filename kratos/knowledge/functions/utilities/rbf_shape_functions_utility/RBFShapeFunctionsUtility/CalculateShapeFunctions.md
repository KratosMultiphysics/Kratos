---
type: C++ Method
title: CalculateShapeFunctions
resource: utilities/rbf_shape_functions_utility.cpp#L38-L83
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/EvaluateRBF
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/sum
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateInverseMultiquadricShapeParameter
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RBFShapeFunctionsUtility::CalculateShapeFunctions( const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, DenseQRPointerType pDenseQR)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [EvaluateRBF](../../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/EvaluateRBF.md)
- [sum](../../../../functions/includes/amatrix_interface/sum.md)
- [CalculateInverseMultiquadricShapeParameter](../../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateInverseMultiquadricShapeParameter.md)