---
type: C++ Function
title: GenerateOuterInclinedRefinedTrimmedBrepSurface
resource: tests/cpp_tests/geometries/test_brep_surface_trimmed.cpp#L390-L461
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve4
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::Pointer GenerateOuterInclinedRefinedTrimmedBrepSurface()`

# Calls

- [GenerateReferenceNurbsSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsSurface.md)
- [KnotsV](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [KnotsU](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [GenerateReferenceNurbsCurve1](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve1.md)
- [GenerateReferenceNurbsCurve2](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve2.md)
- [GenerateReferenceNurbsCurve3](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve3.md)
- [GenerateReferenceNurbsCurve4](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve4.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/KRATOS_TEST_CASE_IN_SUITE.md)