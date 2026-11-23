---
type: Python Method
title: CheckVtmFile
resource: tests/test_vtu_output.py#L540-L565
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CheckVtmFile(test_class: kratos_unittest.TestCase, vtm_file_name: str, vtu_file_name_list: 'list[str]'):`

# Calls

- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)

# Called by

- [test_PointVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)