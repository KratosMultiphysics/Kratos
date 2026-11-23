---
type: Python Method
title: tearDown
resource: tests/test_restart.py#L33-L36
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def tearDown(self):`

# Calls

- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)
- [DeleteDirectoryIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)