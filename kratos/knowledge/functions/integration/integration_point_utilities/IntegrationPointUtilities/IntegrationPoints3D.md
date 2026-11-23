---
type: C++ Method
title: IntegrationPoints3D
resource: integration/integration_point_utilities.cpp#L206-L241
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepVolumeSbmIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::IntegrationPoints3D( typename IntegrationPointsArrayType::iterator& rIntegrationPointsBegin, SizeType PointsInU, SizeType PointsInV, SizeType PointsInW, double U0, double U1, double V0, double V1, double W0, double W1)`

# Calls

- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [CreateIntegrationPoints](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateIntegrationPoints.md)
- [CreateBrepVolumeSbmIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepVolumeSbmIntegrationPoints.md)