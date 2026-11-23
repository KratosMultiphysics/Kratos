---
type: C++ Method
title: PreserveDistancesOnDividedElements
resource: processes/parallel_distance_calculation_process.cpp#L388-L403
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/SetDistancesOnDividedElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelDistanceCalculationProcess<TDim>::PreserveDistancesOnDividedElements( Geometry<Node>& rGeometry, const array_1d<double,TDim+1>& rDist)`

# Called by

- [SetDistancesOnDividedElements](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/SetDistancesOnDividedElements.md)