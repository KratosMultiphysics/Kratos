---
type: Python Method
title: test_file_serialization
resource: tests/test_kratos_parameters.py#L997-L1016
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/WriteJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_file_serialization(self):`

# Calls

- [WriteJsonString](../../../../functions/sources/kratos_parameters/Parameters/WriteJsonString.md)
- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)