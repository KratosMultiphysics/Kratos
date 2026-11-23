---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_tetrahedra_3d_10.cpp#L72-L75
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_tetrahedra_3d_10/GenerateReferenceTetrahedra3D10
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_tetrahedra_3d_10/GenerateCurvedTetrahedra3D10
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateDistanceFrom3DGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateReferenceTetrahedra3D10](../../../../../functions/tests/cpp_tests/geometries/test_tetrahedra_3d_10/GenerateReferenceTetrahedra3D10.md)
- [GenerateCurvedTetrahedra3D10](../../../../../functions/tests/cpp_tests/geometries/test_tetrahedra_3d_10/GenerateCurvedTetrahedra3D10.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)
- [CalculateDistanceFrom3DGeometry](../../../../../functions/utilities/geometry_utilities/CalculateDistanceFrom3DGeometry.md)