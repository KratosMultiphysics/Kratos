---
type: C++ Function
title: GenerateTrimmedBrepSurface
resource: tests/cpp_tests/geometries/test_brep_surface_trimmed.cpp#L104-L129
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateNurbsSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateNurbsCurve1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateNurbsCurve2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateNurbsCurve3
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::Pointer GenerateTrimmedBrepSurface()`

# Calls

- [GenerateNurbsSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateNurbsSurface.md)
- [GenerateNurbsCurve1](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateNurbsCurve1.md)
- [GenerateNurbsCurve2](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateNurbsCurve2.md)
- [GenerateNurbsCurve3](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateNurbsCurve3.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/KRATOS_TEST_CASE_IN_SUITE.md)