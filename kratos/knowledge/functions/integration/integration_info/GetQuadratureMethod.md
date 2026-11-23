---
type: C++ Function
title: GetQuadratureMethod
resource: integration/integration_info.h#L105-L108
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/GetIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/integration/test_integration_info/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`QuadratureMethod GetQuadratureMethod(IndexType DimensionIndex) const`

# Called by

- [CreateQuadraturePointGeometries](../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateQuadraturePointGeometries.md)
- [GetIntegrationMethod](../../../functions/integration/integration_info/IntegrationInfo/GetIntegrationMethod.md)
- [CreateIntegrationPoints1D](../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1D.md)
- [CreateIntegrationPoints2D](../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2D.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/integration/test_integration_info/KRATOS_TEST_CASE_IN_SUITE.md)