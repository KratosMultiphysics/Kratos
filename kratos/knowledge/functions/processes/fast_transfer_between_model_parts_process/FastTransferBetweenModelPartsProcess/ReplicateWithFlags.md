---
type: C++ Method
title: ReplicateWithFlags
resource: processes/fast_transfer_between_model_parts_process.cpp#L353-L458
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/ReorderAllIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FastTransferBetweenModelPartsProcess::ReplicateWithFlags()`

# Calls

- [ReorderAllIds](../../../../functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/ReorderAllIds.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [Execute](../../../../functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/Execute.md)