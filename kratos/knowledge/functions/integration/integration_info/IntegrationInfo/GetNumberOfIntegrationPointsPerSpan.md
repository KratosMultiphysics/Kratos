---
type: C++ Method
title: GetNumberOfIntegrationPointsPerSpan
resource: integration/integration_info.cpp#L104-L107
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/integration/integration_info/IntegrationInfo/GetIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/integration/test_integration_info/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepSurfaceSbmIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepVolumeSbmIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IntegrationInfo::SizeType IntegrationInfo::GetNumberOfIntegrationPointsPerSpan(IndexType DimensionIndex) const`

# Called by

- [GetIntegrationMethod](../../../../functions/integration/integration_info/IntegrationInfo/GetIntegrationMethod.md)
- [CreateIntegrationPoints1D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1D.md)
- [CreateIntegrationPoints2D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2D.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/integration/test_integration_info/KRATOS_TEST_CASE_IN_SUITE.md)
- [CreateBrepSurfaceSbmIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepSurfaceSbmIntegrationPoints.md)
- [CreateBrepVolumeSbmIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepVolumeSbmIntegrationPoints.md)
- [CreateBrepSurfaceTrimmingIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints.md)