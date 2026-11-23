---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_coupling_geometry_quadrature.cpp#L160-L175
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/GenerateNurbsCurveCouplingGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/GenerateNurbsCurveOnSurfaceCouplingGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/SetNumberOfIntegrationPointsPerSpan
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryCurvesCreateIntegrationPoints, KratosCoreCouplingGeometriesFastSuite)`

# Calls

- [GenerateNurbsCurveCouplingGeometry](../../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/GenerateNurbsCurveCouplingGeometry.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [GenerateNurbsCurveOnSurfaceCouplingGeometry](../../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/GenerateNurbsCurveOnSurfaceCouplingGeometry.md)
- [GetNumberOfIntegrationPointsPerSpan](../../../../../functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan.md)
- [SetNumberOfIntegrationPointsPerSpan](../../../../../functions/integration/integration_info/IntegrationInfo/SetNumberOfIntegrationPointsPerSpan.md)