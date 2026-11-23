---
type: C++ Method
title: ComputeIntersectionNormal
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L531-L548
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CorrectDistanceOrientation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeIntersectionNormal( const Element::GeometryType& rGeometry, const Vector& rElementalDistances, array_1d<double,3>& rNormal)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CorrectDistanceOrientation](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CorrectDistanceOrientation.md)