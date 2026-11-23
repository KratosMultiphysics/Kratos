---
type: C++ Method
title: OptimizeOrdering
resource: processes/reorder_and_optimize_modelpart_process.cpp#L170-L211
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/ReorderElements
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReorderAndOptimizeModelPartProcess::OptimizeOrdering()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ReorderElements](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/ReorderElements.md)

# Called by

- [Execute](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/Execute.md)