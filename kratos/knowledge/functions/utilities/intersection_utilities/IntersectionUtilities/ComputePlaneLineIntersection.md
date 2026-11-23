---
type: C++ Method
title: ComputePlaneLineIntersection
resource: utilities/intersection_utilities.cpp#L321-L357
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedGeometryIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int IntersectionUtilities::ComputePlaneLineIntersection( const array_1d<double,3>& rPlaneBasePoint, const array_1d<double,3>& rPlaneNormal, const array_1d<double,3>& rLinePoint1, const array_1d<double,3>& rLinePoint2, array_1d<double,3>& rIntersectionPoint, const double Epsilon)`

# Called by

- [ComputeExtrapolatedGeometryIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedGeometryIntersections.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE.md)