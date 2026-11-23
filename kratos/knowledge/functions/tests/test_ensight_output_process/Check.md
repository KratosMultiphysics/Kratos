---
type: Python Function
title: Check
resource: tests/test_ensight_output_process.py#L173-L195
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def Check(output_file,reference_file, file_format, extension):`

# Calls

- [SetString](../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [AddEmptyValue](../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)