---
type: Python Function
title: ComparePoints
resource: python_scripts/compare_two_files_check_process.py#L418-L433
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareVtkFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ComparePoints(lines_ref, lines_out, line_counter):`

# Calls

- [__CheckCloseValues](../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues.md)

# Called by

- [__CompareVtkFile](../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareVtkFile.md)