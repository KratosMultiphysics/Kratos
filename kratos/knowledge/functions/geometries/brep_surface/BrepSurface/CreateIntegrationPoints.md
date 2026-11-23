---
type: C++ Method
title: CreateIntegrationPoints
resource: geometries/brep_surface.h#L472-L514
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepSurfaceSbmIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateIntegrationPoints( IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo) const override`

# Calls

- [CreateBrepSurfaceSbmIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepSurfaceSbmIntegrationPoints.md)
- [CreateBrepSurfaceTrimmingIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints.md)