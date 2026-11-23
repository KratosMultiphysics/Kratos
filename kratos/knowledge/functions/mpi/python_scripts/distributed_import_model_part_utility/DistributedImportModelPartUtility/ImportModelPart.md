---
type: Python Method
title: ImportModelPart
resource: mpi/python_scripts/distributed_import_model_part_utility.py#L27-L137
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddMissingParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/reorder_consecutive_model_part_io/ReorderConsecutiveModelPartIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/LoadRestart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ImportModelPart(self):`

# Calls

- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [AddMissingParameters](../../../../../functions/sources/kratos_parameters/Parameters/AddMissingParameters.md)
- [GetBool](../../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [ReorderConsecutiveModelPartIO](../../../../../functions/includes/reorder_consecutive_model_part_io/ReorderConsecutiveModelPartIO.md)
- [SetString](../../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [RemoveValue](../../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [AddEmptyValue](../../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetBool](../../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [LoadRestart](../../../../../functions/python_scripts/restart_utility/RestartUtility/LoadRestart.md)