---
type: C++ Method
title: ComputeNegativeSideShapeFunctionsAndWeights
resource: modified_shape_functions/modified_shape_functions.cpp#L110-L132
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeSubdivisions
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndWeights( Matrix &rNegativeSideShapeFunctionsValues, Vector &rNegativeSideWeightsValues, const IntegrationMethodType IntegrationMethod)`

# Calls

- [GetNegativeSubdivisions](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeSubdivisions.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)