---
type: C++ Method
title: CalculateInverseMultiquadricShapeParameter
resource: utilities/rbf_shape_functions_utility.cpp#L146-L167
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(const Matrix& rPoints)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CalculateShapeFunctions](../../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctions.md)