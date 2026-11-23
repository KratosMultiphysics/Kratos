---
type: C++ Method
title: EvaluateRBF
resource: utilities/rbf_shape_functions_utility.cpp#L25-L36
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctionsAndInterpolation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double RBFShapeFunctionsUtility::EvaluateRBF( const double x, const double h)`

# Called by

- [CalculateShapeFunctions](../../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctions.md)
- [CalculateShapeFunctionsAndInterpolation](../../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctionsAndInterpolation.md)