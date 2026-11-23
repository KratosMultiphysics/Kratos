---
type: Python Method
title: __CompareSolMetricFile
resource: python_scripts/compare_two_files_check_process.py#L352-L396
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/ConvertStringToListFloat
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CompareSolMetricFile(self):`

# Calls

- [__GetFileLines](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines.md)
- [ConvertStringToListFloat](../../../../functions/python_scripts/compare_two_files_check_process/ConvertStringToListFloat.md)

# Called by

- [ExecuteFinalize](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize.md)