---
type: Python Method
title: test_nodal_print
resource: tests/test_print_info_in_file.py#L17-L93
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_nodal_print(self):`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [remove](../../../../functions/sources/kratos_filesystem/remove.md)