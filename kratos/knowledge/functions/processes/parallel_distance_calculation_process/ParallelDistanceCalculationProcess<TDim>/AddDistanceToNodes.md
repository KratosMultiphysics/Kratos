---
type: C++ Method
title: AddDistanceToNodes
resource: processes/parallel_distance_calculation_process.cpp#L162-L240
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/ExtendDistancesByLayer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelDistanceCalculationProcess<TDim>::AddDistanceToNodes( Geometry<NodeType>& rGeometry, const BoundedMatrix<double,TDim+1,TDim>& rDN_DX, const double& Volume)`

# Called by

- [ExtendDistancesByLayer](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/ExtendDistancesByLayer.md)