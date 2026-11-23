---
type: Python Method
title: __CompareVtkFile
resource: python_scripts/compare_two_files_check_process.py#L398-L542
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CheckHeader
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/ComparePoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareCells
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareCellTypes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CompareVtkFile(self):`

# Calls

- [__GetFileLines](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__GetFileLines.md)
- [CheckHeader](../../../../functions/python_scripts/compare_two_files_check_process/CheckHeader.md)
- [ComparePoints](../../../../functions/python_scripts/compare_two_files_check_process/ComparePoints.md)
- [CompareCells](../../../../functions/python_scripts/compare_two_files_check_process/CompareCells.md)
- [CompareCellTypes](../../../../functions/python_scripts/compare_two_files_check_process/CompareCellTypes.md)
- [CompareData](../../../../functions/python_scripts/compare_two_files_check_process/CompareData.md)

# Called by

- [ExecuteFinalize](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize.md)