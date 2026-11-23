---
type: C++ Function
title: CalculateExactDistancesToPlane
resource: utilities/geometry_utilities.h#L515-L546
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/CalculateTetrahedraIntersectionPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/parallel_distance_calculation_process/if
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void CalculateExactDistancesToPlane( const GeometryType& rThisGeometry, array_1d<double, TSize>& rDistances )`

# Calls

- [CalculateTetrahedraIntersectionPoints](../../../functions/utilities/geometry_utilities/CalculateTetrahedraIntersectionPoints.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [if](../../../functions/processes/parallel_distance_calculation_process/if.md)
- [Execute](../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)