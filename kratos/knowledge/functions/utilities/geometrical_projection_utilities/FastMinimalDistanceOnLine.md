---
type: C++ Function
title: FastMinimalDistanceOnLine
resource: utilities/geometrical_projection_utilities.h#L235-L251
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceByBruteForce
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/TestFastMinimalDistanceOnLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double FastMinimalDistanceOnLine( const TGeometryType& rGeometry, const Point& rPoint, const double Tolerance = 1.0e-9 )`

# Called by

- [CalculateDistance](../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistance.md)
- [CalculateDistanceByBruteForce](../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceByBruteForce.md)
- [TestFastMinimalDistanceOnLine](../../../functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/TestFastMinimalDistanceOnLine.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/KRATOS_TEST_CASE_IN_SUITE.md)