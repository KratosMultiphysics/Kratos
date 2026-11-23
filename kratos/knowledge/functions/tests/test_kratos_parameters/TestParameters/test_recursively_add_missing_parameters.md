---
type: Python Method
title: test_recursively_add_missing_parameters
resource: tests/test_kratos_parameters.py#L431-L439
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyAddMissingParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_recursively_add_missing_parameters(self): # only missing parameters are added, no complaints if there already exist more than in the defaults`

# Calls

- [RecursivelyAddMissingParameters](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyAddMissingParameters.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)