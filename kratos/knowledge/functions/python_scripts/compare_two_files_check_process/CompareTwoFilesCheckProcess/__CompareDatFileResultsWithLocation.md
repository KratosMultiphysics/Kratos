---
type: Python Method
title: __CompareDatFileResultsWithLocation
resource: python_scripts/compare_two_files_check_process.py#L307-L321
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileVariablesTimeHistory
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CompareDatFileResultsWithLocation(self, lines_ref, lines_out, variable_names):`

# Calls

- [__CheckCloseValues](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues.md)

# Called by

- [__CompareDatFileVariablesTimeHistory](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileVariablesTimeHistory.md)