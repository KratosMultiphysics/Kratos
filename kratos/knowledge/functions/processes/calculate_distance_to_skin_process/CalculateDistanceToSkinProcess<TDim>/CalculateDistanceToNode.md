---
type: C++ Method
title: CalculateDistanceToNode
resource: processes/calculate_distance_to_skin_process.cpp#L201-L236
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/plane_3d/Plane3D/CalculateSignedDistance
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateElementalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CalculateDistanceToSkinProcess<TDim>::CalculateDistanceToNode( Node &rNode, PointerVector<GeometricalObject>& rIntersectedObjects, const double Epsilon)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [CalculateSignedDistance](../../../../functions/geometries/plane_3d/Plane3D/CalculateSignedDistance.md)

# Called by

- [CalculateElementalDistances](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateElementalDistances.md)