---
type: Python Method
title: __get_gidpost_flag
resource: python_scripts/gid_output_process.py#L369-L379
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_InitializeGiDIO
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __get_gidpost_flag(self, param, label, dictionary):`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [_InitializeGiDIO](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_InitializeGiDIO.md)