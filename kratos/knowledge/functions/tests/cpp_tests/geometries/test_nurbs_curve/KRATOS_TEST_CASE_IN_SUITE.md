---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_nurbs_curve.cpp#L168-L205
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve/GenerateReferenceCurve2d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/NumberOfKnots
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve/GenerateCircle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve/GenerateReferenceCurve3d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/SetQuadratureMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/SetNumberOfIntegrationPointsPerSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(NurbsCurve2d, KratosCoreNurbsGeometriesFastSuite)`

# Calls

- [GenerateReferenceCurve2d](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve/GenerateReferenceCurve2d.md)
- [NumberOfKnots](../../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/NumberOfKnots.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [GetT0](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0.md)
- [GetT1](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1.md)
- [GenerateCircle](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve/GenerateCircle.md)
- [GenerateReferenceCurve3d](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve/GenerateReferenceCurve3d.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [SetQuadratureMethod](../../../../../functions/integration/integration_info/IntegrationInfo/SetQuadratureMethod.md)
- [SetNumberOfIntegrationPointsPerSpan](../../../../../functions/integration/integration_info/IntegrationInfo/SetNumberOfIntegrationPointsPerSpan.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)