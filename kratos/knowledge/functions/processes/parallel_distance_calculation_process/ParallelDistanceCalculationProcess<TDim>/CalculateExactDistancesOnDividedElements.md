---
type: C++ Method
title: CalculateExactDistancesOnDividedElements
resource: processes/parallel_distance_calculation_process.cpp#L360-L370
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/CalculateTetrahedraDistances
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/SetDistancesOnDividedElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelDistanceCalculationProcess<TDim>::CalculateExactDistancesOnDividedElements( Geometry<Node>& rGeometry, array_1d<double,TDim+1>& rDist)`

# Calls

- [CalculateTetrahedraDistances](../../../../functions/utilities/geometry_utilities/CalculateTetrahedraDistances.md)

# Called by

- [SetDistancesOnDividedElements](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/SetDistancesOnDividedElements.md)