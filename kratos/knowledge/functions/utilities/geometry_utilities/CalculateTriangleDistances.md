---
type: C++ Function
title: CalculateTriangleDistances
resource: utilities/geometry_utilities.h#L480-L506
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/CalculateTetrahedraIntersectionPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void CalculateTriangleDistances( const GeometryType& rGeometry, array_1d<double, TSize>& rDistances )`

# Calls

- [CalculateTetrahedraIntersectionPoints](../../../functions/utilities/geometry_utilities/CalculateTetrahedraIntersectionPoints.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [PointDistanceToLineSegment3D](../../../functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D.md)

# Called by

- [Execute](../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE.md)