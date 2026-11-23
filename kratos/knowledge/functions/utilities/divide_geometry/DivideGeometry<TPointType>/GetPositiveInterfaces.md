---
type: C++ Method
title: GetPositiveInterfaces
resource: utilities/divide_geometry.cpp#L142-L145
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_divide_tetrahedra_3d_4/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_divide_triangle_2d_3/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_divide_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<typename DivideGeometry<TPointType>::IndexedPointGeometryPointerType> DivideGeometry<TPointType>::GetPositiveInterfaces() const`

# Called by

- [ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues.md)
- [ComputePositiveSideInterfaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_tetrahedra_3d_4/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_triangle_2d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [ComputeElementSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin.md)