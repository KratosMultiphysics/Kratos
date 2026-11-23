---
type: C++ Method
title: GetNegativeSubdivisions
resource: utilities/divide_geometry.cpp#L136-L139
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideDomainSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideShapeFunctionsAndWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeExteriorFaceAreaNormals
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
---

# Signature

`std::vector<typename DivideGeometry<TPointType>::IndexedPointGeometryPointerType> DivideGeometry<TPointType>::GetNegativeSubdivisions() const`

# Called by

- [ComputeNegativeSideDomainSize](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideDomainSize.md)
- [ComputeNegativeSideShapeFunctionsAndWeights](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideShapeFunctionsAndWeights.md)
- [ComputeNegativeSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeSideShapeFunctionsAndGradientsValues.md)
- [ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues.md)
- [ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues.md)
- [ComputeNegativeExteriorFaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeNegativeExteriorFaceAreaNormals.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_tetrahedra_3d_4/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_triangle_2d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE.md)