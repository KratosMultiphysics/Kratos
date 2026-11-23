---
type: C++ Method
title: CalculateShapeFunctionsAndInterpolation
resource: utilities/rbf_shape_functions_utility.cpp#L99-L144
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
---

# Signature

`double RBFShapeFunctionsUtility::CalculateShapeFunctionsAndInterpolation( const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, Vector& rY)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [EvaluateRBF](../../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/EvaluateRBF.md)