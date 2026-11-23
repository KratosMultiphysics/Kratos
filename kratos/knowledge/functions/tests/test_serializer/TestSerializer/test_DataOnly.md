---
type: Python Method
title: test_DataOnly
resource: tests/test_serializer.py#L102-L213
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_serializer/create_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_serializer/comparator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_DataOnly(self):`

# Calls

- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)
- [create_model_part](../../../../functions/tests/test_serializer/create_model_part.md)
- [comparator](../../../../functions/tests/test_serializer/comparator.md)