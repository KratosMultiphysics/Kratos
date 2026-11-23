---
type: C++ Method
title: Execute
resource: processes/reorder_and_optimize_modelpart_process.cpp#L40-L140
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/OptimizeOrdering
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/ActualizeSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReorderAndOptimizeModelPartProcess::Execute()`

# Calls

- [OptimizeOrdering](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/OptimizeOrdering.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetData](../../../../functions/sources/gid_io/GidIOBase/GetData.md)
- [ActualizeSubModelPart](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/ActualizeSubModelPart.md)