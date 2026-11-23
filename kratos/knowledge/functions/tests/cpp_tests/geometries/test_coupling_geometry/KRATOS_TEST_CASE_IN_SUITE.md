---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_coupling_geometry.cpp#L76-L128
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry/GenerateNodeTriangle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry/GenerateNodeTriangle2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/geometrical_object/SetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry/GeneratePointTriangle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry/GeneratePointTriangle3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryNode, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateNodeTriangle](../../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry/GenerateNodeTriangle.md)
- [GenerateNodeTriangle2](../../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry/GenerateNodeTriangle2.md)
- [CreateNewProperties](../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [SetGeometry](../../../../../functions/includes/geometrical_object/SetGeometry.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GeneratePointTriangle](../../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry/GeneratePointTriangle.md)
- [GeneratePointTriangle3D](../../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry/GeneratePointTriangle3D.md)