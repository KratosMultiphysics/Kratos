---
type: C++ Method
title: CreateBrepVolumeSbmIntegrationPoints
resource: utilities/geometry_utilities/brep_sbm_utilities.cpp#L134-L270
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
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints3D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_volume/BrepVolume/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BrepSbmUtilities<TNodeType>::CreateBrepVolumeSbmIntegrationPoints( const std::vector<double>& rSpansU, const std::vector<double>& rSpansV, const std::vector<double>& rSpansW, const GeometrySurrogateArrayType& rSurrogateOuterLoopGeometries, const GeometrySurrogateArrayType& rSurrogateInnerLoopGeometries, IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo)`

# Calls

- [FindKnotSpans1D](../../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/FindKnotSpans1D.md)
- [GetNumberOfIntegrationPointsPerSpan](../../../../../functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan.md)
- [IntegrationPoints3D](../../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints3D.md)

# Called by

- [CreateIntegrationPoints](../../../../../functions/geometries/brep_volume/BrepVolume/CreateIntegrationPoints.md)