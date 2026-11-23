---
type: Python Method
title: __CompareDatFileVariablesTimeHistory
resource: python_scripts/compare_two_files_check_process.py#L258-L283
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
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileResultsWithLocation
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CompareDatFileVariablesTimeHistory(self):`

# Calls

- [__GetFileLines](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines.md)
- [__CompareDatFileComments](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileComments.md)
- [__CompareDatFileResultsWithLocation](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileResultsWithLocation.md)

# Called by

- [ExecuteFinalize](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize.md)