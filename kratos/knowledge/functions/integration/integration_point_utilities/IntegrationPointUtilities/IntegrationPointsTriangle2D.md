---
type: C++ Method
title: IntegrationPointsTriangle2D
resource: integration/integration_point_utilities.cpp#L243-L265
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::IntegrationPointsTriangle2D( typename IntegrationPointsArrayType::iterator& rIntegrationPointsBegin, SizeType PointsIndex, double U0, double U1, double U2, double V0, double V1, double V2)`

# Calls

- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [CreateBrepSurfaceTrimmingIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints.md)