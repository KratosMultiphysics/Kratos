---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_nurbs_volume.cpp#L169-L192
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume/GenerateTruncatedPyramid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume/GenerateDistortedCube
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGeometryIntegrationPoints1, KratosCoreNurbsGeometriesFastSuite)`

# Calls

- [GenerateTruncatedPyramid](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume/GenerateTruncatedPyramid.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [GenerateDistortedCube](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume/GenerateDistortedCube.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)