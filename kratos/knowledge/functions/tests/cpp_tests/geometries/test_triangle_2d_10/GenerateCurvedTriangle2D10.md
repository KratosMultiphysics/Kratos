---
type: C++ Function
title: GenerateCurvedTriangle2D10
resource: tests/cpp_tests/geometries/test_triangle_2d_10.cpp#L54-L65
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_triangle_2d_10/GenerateReferenceNodes2D10
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_triangle_2d_10/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Geometry<NodeType>::Pointer GenerateCurvedTriangle2D10()`

# Calls

- [GenerateReferenceNodes2D10](../../../../../functions/tests/cpp_tests/geometries/test_triangle_2d_10/GenerateReferenceNodes2D10.md)
- [Y0](../../../../../functions/includes/node/Y0.md)
- [X0](../../../../../functions/includes/node/X0.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_triangle_2d_10/KRATOS_TEST_CASE_IN_SUITE.md)