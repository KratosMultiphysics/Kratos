---
type: C++ Function
title: ShapeFunctionSanityCheck
resource: tests/cpp_tests/geometries/test_prism_3d_15.cpp#L204-L276
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_prism_3d_15/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ShapeFunctionSanityCheck( Geometry<Node>::Pointer geom, const GeometryData::IntegrationMethod method, bool isFullIntegration)`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_prism_3d_15/KRATOS_TEST_CASE_IN_SUITE.md)