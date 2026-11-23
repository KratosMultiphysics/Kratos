---
type: C++ Method
title: ComputeDimension
resource: processes/find_conditions_neighbours_process.cpp#L170-L181
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindConditionsNeighboursProcess::ComputeDimension()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [Execute](../../../../functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/Execute.md)