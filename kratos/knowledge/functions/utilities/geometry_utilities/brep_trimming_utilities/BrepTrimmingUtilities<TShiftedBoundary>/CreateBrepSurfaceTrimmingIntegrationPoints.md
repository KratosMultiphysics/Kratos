---
type: C++ Method
title: CreateBrepSurfaceTrimmingIntegrationPoints
resource: utilities/geometry_utilities/brep_trimming_utilities.cpp#L18-L151
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/GetTessellation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/ToIntPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/Triangulate_OPT
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPointsTriangle2D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_surface/BrepSurface/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BrepTrimmingUtilities<TShiftedBoundary>::CreateBrepSurfaceTrimmingIntegrationPoints( IntegrationPointsArrayType& rIntegrationPoints, const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rOuterLoops, const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rInnerLoops, const std::vector<double>& rSpansU, const std::vector<double>& rSpansV, IntegrationInfo& rIntegrationInfo)`

# Calls

- [GetTessellation](../../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/GetTessellation.md)
- [ToIntPoint](../../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/ToIntPoint.md)
- [GetNumberOfIntegrationPointsPerSpan](../../../../../functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan.md)
- [IntegrationPoints2D](../../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D.md)
- [Triangulate_OPT](../../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/Triangulate_OPT.md)
- [IntegrationPointsTriangle2D](../../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPointsTriangle2D.md)

# Called by

- [CreateIntegrationPoints](../../../../../functions/geometries/brep_surface/BrepSurface/CreateIntegrationPoints.md)