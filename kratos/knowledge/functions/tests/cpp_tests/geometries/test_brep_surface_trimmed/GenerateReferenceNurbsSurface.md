---
type: C++ Function
title: GenerateReferenceNurbsSurface
resource: tests/cpp_tests/geometries/test_brep_surface_trimmed.cpp#L131-L159
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateOuterInclinedTrimmedBrepSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateOuterInclinedRefinedTrimmedBrepSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateInnerRectRefinedTrimmedBrepSurface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename NurbsSurfaceGeometry<3, PointerVector<NodeType>>::Pointer GenerateReferenceNurbsSurface()`

# Called by

- [GenerateOuterInclinedTrimmedBrepSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateOuterInclinedTrimmedBrepSurface.md)
- [GenerateOuterInclinedRefinedTrimmedBrepSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateOuterInclinedRefinedTrimmedBrepSurface.md)
- [GenerateInnerRectRefinedTrimmedBrepSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateInnerRectRefinedTrimmedBrepSurface.md)