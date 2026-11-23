---
type: C++ Method
title: ExecuteInitialize
resource: processes/split_internal_interfaces_process.cpp#L20-L69
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/SplitBoundary
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SplitInternalInterfacesProcess::ExecuteInitialize()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [GenericFindElementalNeighboursProcess](../../../../functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess.md)
- [SplitBoundary](../../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/SplitBoundary.md)