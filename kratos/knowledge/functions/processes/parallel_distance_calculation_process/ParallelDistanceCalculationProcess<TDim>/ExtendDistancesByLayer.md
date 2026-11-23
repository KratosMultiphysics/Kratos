---
type: C++ Method
title: ExtendDistancesByLayer
resource: processes/parallel_distance_calculation_process.cpp#L406-L475
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/AddDistanceToNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelDistanceCalculationProcess<TDim>::ExtendDistancesByLayer()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [AddDistanceToNodes](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/AddDistanceToNodes.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)

# Called by

- [Execute](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/Execute.md)