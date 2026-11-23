---
type: C++ Method
title: PointDistanceToLineSegment3D
resource: utilities/geometry_utilities.cpp#L138-L176
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/line_3d_2/Line3D2/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<2>/CalculatePointDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToTwoIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTetrahedraDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTriangleDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double GeometryUtils::PointDistanceToLineSegment3D( const Point& rLinePoint1, const Point& rLinePoint2, const Point& rToPoint )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CalculateDistance](../../../../functions/geometries/line_3d_2/Line3D2/CalculateDistance.md)
- [CalculatePointDistance](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<2>/CalculatePointDistance.md)
- [CalcSignedDistancesToTwoIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToTwoIntNodes.md)
- [CalculateTetrahedraDistances](../../../../functions/utilities/geometry_utilities/CalculateTetrahedraDistances.md)
- [CalculateTriangleDistances](../../../../functions/utilities/geometry_utilities/CalculateTriangleDistances.md)