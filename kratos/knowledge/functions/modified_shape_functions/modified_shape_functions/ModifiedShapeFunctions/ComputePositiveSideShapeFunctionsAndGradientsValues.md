---
type: C++ Method
title: ComputePositiveSideShapeFunctionsAndGradientsValues
resource: modified_shape_functions/modified_shape_functions.cpp#L135-L158
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_ausas_incised_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_ausas_incised_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues( Matrix &rPositiveSideShapeFunctionsValues, ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues, Vector &rPositiveSideWeightsValues, const IntegrationMethodType IntegrationMethod)`

# Calls

- [GetPositiveSubdivisions](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_ausas_incised_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_ausas_incised_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)