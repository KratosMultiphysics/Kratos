---
type: C++ Method
title: CalculateShapeFunctions
resource: utilities/mls_shape_functions_utility.cpp#L148-L209
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/mls_shape_functions_utility/MLSShapeFunctionsUtility/CalculateKernel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quaternion/W-2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MLSShapeFunctionsUtility::CalculateShapeFunctions( const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN)`

# Calls

- [CalculateKernel](../../../../functions/utilities/mls_shape_functions_utility/MLSShapeFunctionsUtility/CalculateKernel.md)
- [W](../../../../functions/utilities/quaternion/W-2.md)