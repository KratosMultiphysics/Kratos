---
type: C++ Method
title: RemoveCommonNodesFromSlaveModelPart
resource: processes/apply_periodic_boundary_condition_process.cpp#L81-L95
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyPeriodicConditionProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyPeriodicConditionProcess::RemoveCommonNodesFromSlaveModelPart()`

# Calls

- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)

# Called by

- [ApplyPeriodicConditionProcess](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyPeriodicConditionProcess.md)