---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_brep_surface_trimmed.cpp#L540-L557
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateTrimmedBrepSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
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

`KRATOS_TEST_CASE_IN_SUITE(NurbsTrimmedBrepSurfaceIntegration, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateTrimmedBrepSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateTrimmedBrepSurface.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [GenerateOuterInclinedTrimmedBrepSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateOuterInclinedTrimmedBrepSurface.md)
- [GenerateOuterInclinedRefinedTrimmedBrepSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateOuterInclinedRefinedTrimmedBrepSurface.md)
- [GenerateInnerRectRefinedTrimmedBrepSurface](../../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateInnerRectRefinedTrimmedBrepSurface.md)