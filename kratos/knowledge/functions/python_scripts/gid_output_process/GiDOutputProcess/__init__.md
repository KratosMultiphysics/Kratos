---
type: Python Method
title: __init__
resource: python_scripts/gid_output_process.py#L94-L133
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_ValidateFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self,model_part,file_name,param = None):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [_ValidateFileName](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_ValidateFileName.md)
- [AddString](../../../../functions/sources/kratos_parameters/Parameters/AddString.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)