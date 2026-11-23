---
type: Python Method
title: tearDown
resource: tests/test_sequential_orchestrator.py#L82-L88
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def tearDown(self): # Remove the checkpointing files`

# Calls

- [DeleteDirectoryIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)