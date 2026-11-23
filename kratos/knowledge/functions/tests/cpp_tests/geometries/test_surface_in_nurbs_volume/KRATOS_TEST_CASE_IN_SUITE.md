---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_surface_in_nurbs_volume.cpp#L139-L214
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_surface_in_nurbs_volume/GenerateCubeNurbsVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_surface_in_nurbs_volume/GenerateCuboidNurbsVolume
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(SurfaceInVolumeGeometryTriangleInCubeTest, KratosCoreNurbsGeometriesFastSuite)`

# Calls

- [GenerateCubeNurbsVolume](../../../../../functions/tests/cpp_tests/geometries/test_surface_in_nurbs_volume/GenerateCubeNurbsVolume.md)
- [GetDefaultIntegrationMethod](../../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [GenerateCuboidNurbsVolume](../../../../../functions/tests/cpp_tests/geometries/test_surface_in_nurbs_volume/GenerateCuboidNurbsVolume.md)