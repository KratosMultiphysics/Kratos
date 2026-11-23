---
type: C++ Function
title: CalculateTetrahedraDistances
resource: utilities/geometry_utilities.h#L433-L472
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
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/CalculateExactDistancesOnDividedElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void CalculateTetrahedraDistances( const GeometryType& rGeometry, array_1d<double, TSize>& rDistances)`

# Calls

- [CalculateTetrahedraIntersectionPoints](../../../functions/utilities/geometry_utilities/CalculateTetrahedraIntersectionPoints.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [PointDistanceToLineSegment3D](../../../functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D.md)

# Called by

- [CalculateExactDistancesOnDividedElements](../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/CalculateExactDistancesOnDividedElements.md)
- [Execute](../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)
- [CalculateDistances3D](../../../functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances3D.md)