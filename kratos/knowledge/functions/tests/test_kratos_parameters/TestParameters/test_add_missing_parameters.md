---
type: Python Method
title: test_add_missing_parameters
resource: tests/test_kratos_parameters.py#L420-L429
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddMissingParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_add_missing_parameters(self): # only missing parameters are added, no complaints if there already exist more than in the defaults`

# Calls

- [AddMissingParameters](../../../../functions/sources/kratos_parameters/Parameters/AddMissingParameters.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)