---
type: Python Method
title: __CompareCSVFile
resource: python_scripts/compare_two_files_check_process.py#L244-L256
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileComments
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDelimittedFileResults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CompareCSVFile(self):`

# Calls

- [__GetFileLines](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines.md)
- [__CompareDatFileComments](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileComments.md)
- [__CompareDelimittedFileResults](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDelimittedFileResults.md)

# Called by

- [ExecuteFinalize](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize.md)