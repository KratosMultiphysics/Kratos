---
type: C++ Method
title: SetLinearStrategy
resource: processes/calculate_embedded_nodal_variable_from_skin_process.h#L610-L633
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetLinearStrategy()`

# Calls

- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)

# Called by

- [Execute](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/Execute.md)