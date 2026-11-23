---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_hexahedra_3d_8.cpp#L104-L111
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_hexahedra_3d_8/GenerateOriginCenterLen1Hexahedra3D8
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_hexahedra_3d_8/GenerateOriginCenterLen2Hexahedra3D8
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_hexahedra_3d_8/GenerateDeformedCenterLen1Hexahedra3D8
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateDistanceFrom3DGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateOriginCenterLen1Hexahedra3D8](../../../../../functions/tests/cpp_tests/geometries/test_hexahedra_3d_8/GenerateOriginCenterLen1Hexahedra3D8.md)
- [GenerateOriginCenterLen2Hexahedra3D8](../../../../../functions/tests/cpp_tests/geometries/test_hexahedra_3d_8/GenerateOriginCenterLen2Hexahedra3D8.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)
- [GenerateDeformedCenterLen1Hexahedra3D8](../../../../../functions/tests/cpp_tests/geometries/test_hexahedra_3d_8/GenerateDeformedCenterLen1Hexahedra3D8.md)
- [CalculateDistanceFrom3DGeometry](../../../../../functions/utilities/geometry_utilities/CalculateDistanceFrom3DGeometry.md)