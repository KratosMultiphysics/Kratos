---
type: C++ Method
title: IntegrationPoints2D
resource: integration/integration_point_utilities.cpp#L174-L204
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2DGauss
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepSurfaceSbmIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::IntegrationPoints2D( typename IntegrationPointsArrayType::iterator& rIntegrationPointsBegin, SizeType PointsInU, SizeType PointsInV, double U0, double U1, double V0, double V1)`

# Calls

- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [CreateIntegrationPoints](../../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CreateIntegrationPoints.md)
- [CreateIntegrationPoints2DGauss](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2DGauss.md)
- [CreateBrepSurfaceSbmIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepSurfaceSbmIntegrationPoints.md)
- [CreateBrepSurfaceTrimmingIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints.md)