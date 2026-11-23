---
type: C++ Method
title: GetPositiveSubdivisions
resource: utilities/divide_geometry.cpp#L130-L133
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/internals/skin_intersection/TrimHigher
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/skin_intersection/TrimLower
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideDomainSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideShapeFunctionsAndWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveExteriorFaceAreaNormals
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

`std::vector<typename DivideGeometry<TPointType>::IndexedPointGeometryPointerType> DivideGeometry<TPointType>::GetPositiveSubdivisions() const`

# Called by

- [TrimHigher](../../../../functions/modeler/internals/skin_intersection/TrimHigher.md)
- [TrimLower](../../../../functions/modeler/internals/skin_intersection/TrimLower.md)
- [ComputePositiveSideDomainSize](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideDomainSize.md)
- [ComputePositiveSideShapeFunctionsAndWeights](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideShapeFunctionsAndWeights.md)
- [ComputePositiveSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideShapeFunctionsAndGradientsValues.md)
- [ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues.md)
- [ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues.md)
- [ComputePositiveExteriorFaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveExteriorFaceAreaNormals.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_tetrahedra_3d_4/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_triangle_2d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE.md)