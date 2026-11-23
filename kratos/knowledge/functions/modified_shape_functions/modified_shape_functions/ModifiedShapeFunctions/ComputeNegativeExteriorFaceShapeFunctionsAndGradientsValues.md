---
type: C++ Method
title: ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues
resource: modified_shape_functions/modified_shape_functions.cpp#L282-L318
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeSubdivisions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeFaceValuesOnOneSide
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE
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

`void ModifiedShapeFunctions::ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues( Matrix &rNegativeExteriorFaceShapeFunctionsValues, ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues, Vector &rNegativeExteriorFaceWeightsValues, const unsigned int FaceId, const IntegrationMethodType IntegrationMethod)`

# Calls

- [GetNegativeSubdivisions](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeSubdivisions.md)
- [ComputeFaceValuesOnOneSide](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeFaceValuesOnOneSide.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)