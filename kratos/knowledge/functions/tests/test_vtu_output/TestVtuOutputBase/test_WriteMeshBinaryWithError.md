---
type: Python Method
title: test_WriteMeshBinaryWithError
resource: tests/test_vtu_output.py#L91-L97
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/WriteVtu
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_WriteMeshBinaryWithError(self):`

# Calls

- [WriteVtu](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/WriteVtu.md)
- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)