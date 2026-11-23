---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_hexahedra_3d_20.cpp#L91-L96
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_hexahedra_3d_20/GenerateCanonicalHexahedra3D20
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateCanonicalHexahedra3D20](../../../../../functions/tests/cpp_tests/geometries/test_hexahedra_3d_20/GenerateCanonicalHexahedra3D20.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)