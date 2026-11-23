---
type: C++ Function
title: ComputeShortestLineBetweenTwoLines
resource: utilities/intersection_utilities.h#L460-L559
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/GeometricalProjectionUtilities/FastMinimalDistanceOnLineWithRadius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static PointerVector<Point> ComputeShortestLineBetweenTwoLines( const TGeometryType& rSegment1, const TGeometryType& rSegment2 )`

# Calls

- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [FastMinimalDistanceOnLineWithRadius](../../../functions/utilities/geometrical_projection_utilities/GeometricalProjectionUtilities/FastMinimalDistanceOnLineWithRadius.md)