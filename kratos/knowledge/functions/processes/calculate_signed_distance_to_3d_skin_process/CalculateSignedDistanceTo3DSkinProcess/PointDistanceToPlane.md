---
type: C++ Method
title: PointDistanceToPlane
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1438-L1454
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToOneIntNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToTwoIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double PointDistanceToPlane( Point& planeBasePoint, array_1d<double, 3>& planeNormal, Point& ToPoint)`

# Called by

- [CalcSignedDistancesToOneIntNode](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToOneIntNode.md)
- [CalcSignedDistancesToTwoIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToTwoIntNodes.md)
- [CalcSignedDistancesToThreeIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToThreeIntNodes.md)
- [CalcSignedDistancesToMoreThanThreeIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes.md)