---
type: C++ Method
title: CreateBrepSurfaceSbmIntegrationPoints
resource: utilities/geometry_utilities/brep_sbm_utilities.cpp#L19-L130
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/FindKnotSpans1D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_surface/BrepSurface/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BrepSbmUtilities<TNodeType>::CreateBrepSurfaceSbmIntegrationPoints( const std::vector<double>& rSpansU, const std::vector<double>& rSpansV, const GeometrySurrogateArrayType& rSurrogateOuterLoopGeometries, const GeometrySurrogateArrayType& rSurrogateInnerLoopGeometries, IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo)`

# Calls

- [FindKnotSpans1D](../../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/FindKnotSpans1D.md)
- [GetNumberOfIntegrationPointsPerSpan](../../../../../functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan.md)
- [IntegrationPoints2D](../../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D.md)

# Called by

- [CreateIntegrationPoints](../../../../../functions/geometries/brep_surface/BrepSurface/CreateIntegrationPoints.md)