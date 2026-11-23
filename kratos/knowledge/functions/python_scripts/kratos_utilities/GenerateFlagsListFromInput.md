---
type: Python Function
title: GenerateFlagsListFromInput
resource: python_scripts/kratos_utilities.py#L118-L121
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GenerateFlagsListFromInput(param):`

# Calls

- [GetStringArray](../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)

# Called by

- [ExecuteInitialize](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize.md)