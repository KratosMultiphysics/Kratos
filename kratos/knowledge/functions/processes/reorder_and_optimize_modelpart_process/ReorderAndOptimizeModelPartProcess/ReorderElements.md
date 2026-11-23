---
type: C++ Method
title: ReorderElements
resource: processes/reorder_and_optimize_modelpart_process.cpp#L213-L236
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/OptimizeOrdering
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReorderAndOptimizeModelPartProcess::ReorderElements()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [OptimizeOrdering](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/OptimizeOrdering.md)