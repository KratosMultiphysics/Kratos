---
type: Python Method
title: test_model_serialization
resource: tests/test_model.py#L77-L98
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_model/TestModel/_create_and_save_model
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_model_serialization(self):`

# Calls

- [_create_and_save_model](../../../../functions/tests/test_model/TestModel/_create_and_save_model.md)
- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)
- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)