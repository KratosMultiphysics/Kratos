---
type: Python Method
title: test_sequential_orchestrator
resource: tests/test_sequential_orchestrator.py#L59-L80
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/exists
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_sequential_orchestrator(self): # Parse sequential orchestrator settings`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [exists](../../../../functions/sources/kratos_filesystem/exists.md)