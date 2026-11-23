---
type: C++ Method
title: ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues
resource: modified_shape_functions/modified_shape_functions.cpp#L215-L240
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeFaceValuesOnOneSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeInterfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeSubdivisions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeInterfacesParentIds
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
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GetModifiedShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues( Matrix &rInterfaceNegativeSideShapeFunctionsValues, ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues, Vector &rInterfaceNegativeSideWeightsValues, const IntegrationMethodType IntegrationMethod)`

# Calls

- [ComputeFaceValuesOnOneSide](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeFaceValuesOnOneSide.md)
- [GetNegativeInterfaces](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeInterfaces.md)
- [GetNegativeSubdivisions](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeSubdivisions.md)
- [GetNegativeInterfacesParentIds](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeInterfacesParentIds.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_ausas_incised_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_tetrahedra_3d_4_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_ausas_incised_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_ausas_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)
- [GetModifiedShapeFunctionsValues](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GetModifiedShapeFunctionsValues.md)