---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_prism_3d_15.cpp#L57-L61
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_prism_3d_15/GenerateRegularPrism3D15
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_prism_3d_15/ShapeFunctionSanityCheck
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Prism3D15EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateRegularPrism3D15](../../../../../functions/tests/cpp_tests/geometries/test_prism_3d_15/GenerateRegularPrism3D15.md)
- [CalculateAreaByIntegration](../../../../../functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration.md)
- [ShapeFunctionSanityCheck](../../../../../functions/tests/cpp_tests/geometries/test_prism_3d_15/ShapeFunctionSanityCheck.md)