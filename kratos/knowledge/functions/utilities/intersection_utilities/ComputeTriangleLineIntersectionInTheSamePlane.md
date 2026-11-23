---
type: C++ Function
title: ComputeTriangleLineIntersectionInTheSamePlane
resource: utilities/intersection_utilities.h#L203-L324
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/intersection_utilities/ComputeTetrahedraLineIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static int ComputeTriangleLineIntersectionInTheSamePlane( //static IntersectionUtilitiesTetrahedraLineIntersectionStatus ComputeTriangleLineIntersectionInTheSamePlane( const TGeometryType& rTriangleGeometry, const TCoordinatesType& rLinePoint1, const TCoordinatesType& rLinePoint2, TCoordinatesType& rIntersectionPoint1, TCoordinatesType& rIntersectionPoint2, int& rSolution,//IntersectionUtilitiesTetrahedraLineIntersectionStatus& rSolution, const double Epsilon = 1e-12 )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [ComputeTetrahedraLineIntersection](../../../functions/utilities/intersection_utilities/ComputeTetrahedraLineIntersection.md)