---
type: C++ Method
title: Execute
resource: processes/generic_find_elements_neighbours_process.cpp#L23-L61
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess/CheckForNeighbourElems
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenericFindElementalNeighboursProcess::Execute()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CheckForNeighbourElems](../../../../functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess/CheckForNeighbourElems.md)