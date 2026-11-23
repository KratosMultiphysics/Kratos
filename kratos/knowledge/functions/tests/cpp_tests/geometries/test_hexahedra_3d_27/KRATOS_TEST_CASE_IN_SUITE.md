---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_hexahedra_3d_27.cpp#L68-L73
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_hexahedra_3d_27/GenerateCanonicalHexahedra3D27
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/Det3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateCanonicalHexahedra3D27](../../../../../functions/tests/cpp_tests/geometries/test_hexahedra_3d_27/GenerateCanonicalHexahedra3D27.md)
- [Det3](../../../../../functions/utilities/math_utils/Det3.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)