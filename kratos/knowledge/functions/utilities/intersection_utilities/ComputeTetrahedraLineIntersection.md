---
type: C++ Function
title: ComputeTetrahedraLineIntersection
resource: utilities/intersection_utilities.h#L362-L446
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/ComputeTriangleLineIntersectionInTheSamePlane
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static int ComputeTetrahedraLineIntersection(//static IntersectionUtilitiesTetrahedraLineIntersectionStatus ComputeTetrahedraLineIntersection( const TGeometryType& rTetrahedraGeometry, const TCoordinatesType& rLinePoint1, const TCoordinatesType& rLinePoint2, TCoordinatesType& rIntersectionPoint1, TCoordinatesType& rIntersectionPoint2, const double Epsilon = 1e-12 )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [ComputeTriangleLineIntersectionInTheSamePlane](../../../functions/utilities/intersection_utilities/ComputeTriangleLineIntersectionInTheSamePlane.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE.md)