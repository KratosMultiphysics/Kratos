---
type: C++ Method
title: Execute
resource: processes/find_conditions_neighbours_process.cpp#L59-L111
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/ComputeDimension
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/CheckForNeighbourFaces
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindConditionsNeighboursProcess::Execute()`

# Calls

- [ComputeDimension](../../../../functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/ComputeDimension.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CheckForNeighbourFaces](../../../../functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/CheckForNeighbourFaces.md)