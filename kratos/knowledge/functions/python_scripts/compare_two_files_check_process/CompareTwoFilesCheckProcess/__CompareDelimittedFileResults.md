---
type: Python Method
title: __CompareDelimittedFileResults
resource: python_scripts/compare_two_files_check_process.py#L299-L305
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareCSVFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CompareDelimittedFileResults(self, lines_ref, lines_out, delimiter):`

# Calls

- [__CheckCloseValues](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues.md)

# Called by

- [__CompareDatFile](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFile.md)
- [__CompareCSVFile](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareCSVFile.md)