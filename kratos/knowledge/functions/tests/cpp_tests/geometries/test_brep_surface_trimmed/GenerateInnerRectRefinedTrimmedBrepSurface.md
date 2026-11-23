---
type: C++ Function
title: GenerateInnerRectRefinedTrimmedBrepSurface
resource: tests/cpp_tests/geometries/test_brep_surface_trimmed.cpp#L463-L537
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
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve5
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve6
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve7
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve8
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve9
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve10
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve11
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve12
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::Pointer GenerateInnerRectRefinedTrimmedBrepSurface()`

# Calls

- [GenerateReferenceNurbsSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsSurface.md)
- [KnotsV](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [KnotsU](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [GenerateReferenceNurbsCurve5](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve5.md)
- [GenerateReferenceNurbsCurve6](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve6.md)
- [GenerateReferenceNurbsCurve7](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve7.md)
- [GenerateReferenceNurbsCurve8](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve8.md)
- [GenerateReferenceNurbsCurve9](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve9.md)
- [GenerateReferenceNurbsCurve10](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve10.md)
- [GenerateReferenceNurbsCurve11](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve11.md)
- [GenerateReferenceNurbsCurve12](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateReferenceNurbsCurve12.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/KRATOS_TEST_CASE_IN_SUITE.md)