---
type: Python Method
title: __CompareCaseFile
resource: python_scripts/compare_two_files_check_process.py#L544-L696
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/ParseCaseFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CompareCaseFile(self):`

# Calls

- [__GetFileLines](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines.md)
- [ParseCaseFile](../../../../functions/python_scripts/compare_two_files_check_process/ParseCaseFile.md)
- [__CheckCloseValues](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues.md)

# Called by

- [ExecuteFinalize](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize.md)