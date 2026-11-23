---
type: Python Method
title: __ComparePostResFileResultsBlock
resource: python_scripts/compare_two_files_check_process.py#L201-L228
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__ComparePostResFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __ComparePostResFileResultsBlock(self, lines1, lines2, current_index):`

# Calls

- [__CheckCloseValues](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues.md)

# Called by

- [__ComparePostResFile](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__ComparePostResFile.md)